import scanpy as sc
import pandas as pd
import pickle
from joblib import Parallel, delayed
import numpy as np
from tqdm import tqdm
import warnings
import scipy.sparse as sp

from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from statsmodels.stats.multitest import multipletests

utils = importr('utils')
base = importr('base')
splines = importr('splines')
stats = importr('stats')
lmtest = importr('lmtest')
fitdistrplus = importr('fitdistrplus')
dplyr = importr('dplyr')
lme4 = importr('lme4')

def _as_dense_vector(x):
    """Return a dense 1D numpy array from dense or sparse AnnData slices."""
    if sp.issparse(x):
        x = x.toarray()
    return np.asarray(x).ravel()


def _fallback_dse_result(n_cells):
    """Fail-safe DSE result for a single gene/dimension fit failure."""
    return np.array([-1.0], dtype=np.float64), np.zeros(n_cells, dtype=np.float64), np.zeros(n_cells, dtype=np.float64)


def _warn_if_k_is_large(adata_list, k, multi=False, cell_type=None, group_id=None):
    cell_counts = [adata.shape[0] for adata in adata_list]
    min_cells = min(cell_counts) if cell_counts else 0
    if min_cells and k > max(10, min_cells // 10):
        warnings.warn(
            f"SpaceExpress_DSE k={k} may be large for the smallest sample ({min_cells} cells). "
            "For small or multi-replicate datasets, consider lowering k to reduce rank-deficient spline fits.",
            RuntimeWarning,
            stacklevel=2,
        )
    if multi and group_id is not None:
        groups = sorted(set(group_id))
        if len(groups) != 2:
            raise ValueError("Current multi-replicate DSE implementation supports exactly two groups.")
        for group in groups:
            if sum(1 for g in group_id if g == group) < 1:
                raise ValueError(f"Group {group!r} has no replicate.")
    if cell_type is not None:
        missing = [i for i, adata in enumerate(adata_list) if cell_type not in adata.obs]
        if missing:
            raise ValueError(f"cell_type column {cell_type!r} is missing from adata_list indices {missing}.")

def spline(df, k):
    """
    Fitting a spline model to the data and calculating the likelihood ratio test statistic
    
    Parameters:
    df (pd.DataFrame): A pandas DataFrame containing the data
    knots (int): The number of knots to use in the spline model
    
    Returns:
    test_statistics (np.array): An array of likelihood ratio test statistics
    """
    # Convert pandas DataFrame to R DataFrame
    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_df = robjects.conversion.py2rpy(df)

    # Assign the DataFrame to an R variable
    robjects.globalenv['input'] = r_df
    robjects.globalenv['df'] = k


    r_script = """
    suppressPackageStartupMessages(library(splines))
    suppressPackageStartupMessages(library(lmtest))
    suppressPackageStartupMessages(library(dplyr))
    
    colnames(input) = c("embedding","gene","group")
    data = input
    data$group = as.factor(data$group)

    rm_index = which(data$gene >= mean(data$gene) + 4*sd(data$gene))
    
    if (length(rm_index)!= 0) {
        data = data[-rm_index,] # removing cells with over 4 sd away from mean 
    }
    zero_variation_groups <- data %>% group_by(group) %>% summarize(variation = var(gene)) %>% filter(variation == 0)

    if (nrow(zero_variation_groups) > 0) {
        chisq = -1
        predictions = rep(0, nrow(input))
        interaction = rep(0, nrow(input))
    } else {

        data_temp = data %>% group_by(group) %>% reframe(gene_scaled = scale(gene,center = F,scale = sd(gene))) %>% ungroup()
        data$gene = data_temp$gene_scaled  # scaling but not centering expression within each replicate

        spline_df = ns(data$embedding, df = df)
        data = data.frame(spline_df, group = as.factor(data$group), gene = data$gene)

        # Running the linear models
        ## model with interaction terms
        formula = as.formula(paste0("gene~ (",paste0("X",seq(1,df),collapse = " + "), ")* group "))
        mat = model.matrix(object = formula, data = data)
        zero_index = which(colSums(mat) == 0)
        if (length(zero_index) > 0) {
          mat_new = mat[,-zero_index]
        }else{
          mat_new = mat
        }

        data = data.frame(gene = data$gene, mat_new)
        full = lm(formula = gene~.+0,data = data)
        predicted_full = as.vector(predict(full))
        while (sum(is.na(coef(full))) > 0) {
            # Get the names of the coefficients
            coeff_names <- names(coef(full))

            # Remove the coefficients with NA values from the formula
            valid_coeffs <- which(!is.na(coef(full)))
            mat_new = mat_new[,valid_coeffs]
            data = data.frame(gene = data$gene, mat_new)
            # Refit the model without the NA coefficients
            full <- lm(formula = gene~.+0, data = data)
        }


        # Check again for NA coefficients after refitting
        if (sum(is.na(summary(full)[["coefficients"]][, 1])) > 0) {
            chisq = -1
            predictions = rep(0, nrow(input))
            interaction = rep(0, nrow(input))

        } else {
            predicted_full = as.vector(predict(full))

            # evaluating the interaction terms
            interaction_index = grep(pattern = ".group",x = colnames(model.matrix(full)))
            interaction_coef = summary(full)[["coefficients"]][interaction_index,1]
            present_main = unlist(strsplit(names(interaction_coef),split = ".group1"))
            interaction_df = model.matrix(full)[,present_main]
            interaction_estimates = interaction_df%*%interaction_coef + summary(full)[["coefficients"]][which(rownames(summary(full)[["coefficients"]]) == "group1"),1]


            ## model without interaction terms
            # formula1 = as.formula(paste0("gene~ ",paste0("X",seq(1,df),collapse = " + "), " + group"))
            mat_null = mat_new[,-interaction_index]
            data = data.frame(gene = data$gene, mat_null)
            null = lm(gene~.+0, data = data)
            ## Likelihood Ratio Test 
            vals <- lrtest(null, full)$Chisq[2]

            # Return the test statistics, predicted values, and interaction estimates
            chisq = vals

            generate_vector <- function(values, indices) {
              length_result <- length(values) + length(indices)
              result <- rep(0, length_result)
              value_positions <- setdiff(seq_along(result), indices)
              result[value_positions] <- values
              return(result)
            }    

            predictions = generate_vector(predicted_full, rm_index)
            interaction = generate_vector(interaction_estimates, rm_index)
        }
    }
    """

    robjects.r(r_script)

    test_statistics = robjects.globalenv['chisq']
    test_statistics = np.array(test_statistics).astype(np.float64)

    predictions = robjects.globalenv['predictions']
    predictions = np.array(predictions).astype(np.float64)

    interaction = robjects.globalenv['interaction']
    interaction = np.array(interaction).astype(np.float64)
        
    return test_statistics, predictions, interaction

def spline_multi_rep(df, k):
    """
    Fitting a spline model to the data and calculating the likelihood ratio test statistic
    
    Parameters:
    df (pd.DataFrame): A pandas DataFrame containing the data
    knots (int): The number of knots to use in the spline model
    
    Returns:
    test_statistics (np.array): An array of likelihood ratio test statistics
    """
    # Convert pandas DataFrame to R DataFrame
    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_df = robjects.conversion.py2rpy(df)

    # Assign the DataFrame to an R variable
    robjects.globalenv['input'] = r_df
    robjects.globalenv['df'] = k

    r_script = """
    suppressPackageStartupMessages(library(splines))
    suppressPackageStartupMessages(library(lmtest))
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(lme4))
    
    colnames(input) = c("embedding","gene","group", "rep")
    data = input
    data$group = as.factor(data$group)
    
    val1 = quantile(data$gene,0.99)
    val2 = mean(data$gene) + 4*sd(data$gene)
    val = min(val1, val2)
    
    rm_index = which(data$gene >= val)
    
    if (length(rm_index)!= 0) {
        data = data[-rm_index,]  # removing cells with over 5 sd away from mean
    }
    
    zero_variation_groups <- data %>% group_by(rep,group) %>% summarize(variation = var(gene), .groups = 'drop') %>% filter(variation == 0)
    
    if (nrow(zero_variation_groups) > 0) {
        chisq = -1
        predictions = rep(0, nrow(input))
        interaction = rep(0, nrow(input))
    } else {
        data_temp = data %>% group_by(rep) %>% reframe(gene_scaled = scale(gene,center = F,scale = sd(gene))) %>% ungroup()
        data$gene = data_temp$gene_scaled  # scaling but not centering expression within each replicate
        spline_df = ns(data$embedding, df = df)
        data = data.frame(spline_df, group = as.factor(data$group), gene = data$gene, rep = data$rep)
       
        model_result <- tryCatch({
            # Running the linear models with modified control to avoid singular fit warnings
            formula = as.formula(paste0("gene ~ (", paste0("X", seq(1, df), collapse = " + "), ") * group  + (1 | rep)"))
            full = lmer(formula = formula, data = data, control = lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE))

            # Running the null model
            formula1 = as.formula(paste0("gene ~ ", paste0("X", seq(1, df), collapse = " + "), " + group  + (1 | rep)"))
            null = lmer(formula = formula1, data = data, control = lmerControl(check.conv.singular = "ignore"))

            X_full = model.matrix(full)
            beta_full = fixef(full)
            full_terms = intersect(colnames(X_full), names(beta_full))
            predicted_full = as.vector(X_full[, full_terms, drop = FALSE] %*% beta_full[full_terms])

            # Evaluating the interaction terms by coefficient name so rank-dropped columns are ignored safely.
            interaction_terms = grep(pattern = ":", x = names(beta_full), value = TRUE)
            interaction_terms = interaction_terms[grepl("group", interaction_terms)]
            if (length(interaction_terms) > 0) {
                interaction_df = X_full[, interaction_terms, drop = FALSE]
                interaction_coef = beta_full[interaction_terms]
                interaction_estimates_all = as.vector(interaction_df %*% interaction_coef)
            } else {
                interaction_estimates_all = rep(0, nrow(data))
            }
            group_terms = grep(pattern = "^group", x = names(beta_full), value = TRUE)
            if (length(group_terms) > 0) {
                interaction_estimates_all = interaction_estimates_all + beta_full[group_terms[1]]
            }

            # Likelihood Ratio Test
            vals <- lrtest(null, full)$Chisq[2]
            list(chisq = vals, predicted_full = predicted_full, interaction_estimates_all = interaction_estimates_all)
        }, error = function(e) {
            list(chisq = -1, predicted_full = rep(0, nrow(data)), interaction_estimates_all = rep(0, nrow(data)))
        })

        chisq = model_result$chisq
        generate_vector <- function(values, indices) {
            length_result <- length(values) + length(indices)
            result <- rep(0, length_result)
            value_positions <- setdiff(seq_along(result), indices)
            result[value_positions] <- values
            return(result)
        }    
        predictions = generate_vector(model_result$predicted_full, rm_index)
        interaction = generate_vector(model_result$interaction_estimates_all, rm_index)
    }
    """

    robjects.r(r_script)

    test_statistics = robjects.globalenv['chisq']
    test_statistics = np.array(test_statistics).astype(np.float64)

    predictions = robjects.globalenv['predictions']
    predictions = np.array(predictions).astype(np.float64)

    interaction = robjects.globalenv['interaction']
    interaction = np.array(interaction).astype(np.float64)
        
    return test_statistics, predictions, interaction

def spline_multi_rep_ct(df, k):
    """
    Fitting a spline model to the data and calculating the likelihood ratio test statistic
    
    Parameters:
    df (pd.DataFrame): A pandas DataFrame containing the data
    knots (int): The number of knots to use in the spline model
    
    Returns:
    test_statistics (np.array): An array of likelihood ratio test statistics
    """
    # Convert pandas DataFrame to R DataFrame
    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_df = robjects.conversion.py2rpy(df)

    # Assign the DataFrame to an R variable
    robjects.globalenv['input'] = r_df
    robjects.globalenv['df'] = k

    r_script = """
    suppressPackageStartupMessages(library(splines))
    suppressPackageStartupMessages(library(lmtest))
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(lme4))
    
    colnames(input) = c("embedding","gene","group", "rep","cell_type")
    data = input
    data$group = as.factor(data$group)
    data$cell_type = as.factor(data$cell_type)
    
    val1 = quantile(data$gene,0.99)
    val2 = mean(data$gene) + 4*sd(data$gene)
    val = min(val1, val2)
    
    rm_index = which(data$gene >= val)
    
    if (length(rm_index)!= 0) {
        data = data[-rm_index,]  # removing cells with over 5 sd away from mean
    }
    
    zero_variation_groups <- data %>% group_by(rep,group) %>% summarize(variation = var(gene), .groups = 'drop') %>% filter(variation == 0)
    
    if (nrow(zero_variation_groups) > 0) {
        chisq = -1
        predictions = rep(0, nrow(input))
        interaction = rep(0, nrow(input))
    } else {
        data_temp = data %>% group_by(rep) %>% reframe(gene_scaled = scale(gene,center = F,scale = sd(gene))) %>% ungroup()
        data$gene = data_temp$gene_scaled  # scaling but not centering expression within each replicate
        spline_df = ns(data$embedding, df = df)
        data = data.frame(spline_df, group = as.factor(data$group), gene = data$gene, rep = data$rep, cell_type = data$cell_type)

        model_result <- tryCatch({
            # Running the linear models with modified control to avoid singular fit warnings
            formula = as.formula(paste0("gene ~ (", paste0("X", seq(1, df), collapse = " + "), ") * group + cell_type + (1 | rep)"))
            full = lmer(formula = formula, data = data, control = lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE))

            # Running the null model
            formula1 = as.formula(paste0("gene ~ ", paste0("X", seq(1, df), collapse = " + "), " + group + cell_type + (1 | rep)"))
            null = lmer(formula = formula1, data = data, control = lmerControl(check.conv.singular = "ignore"))

            X_full = model.matrix(full)
            beta_full = fixef(full)
            non_cell_terms = setdiff(intersect(colnames(X_full), names(beta_full)), grep(pattern = "^cell_type", x = names(beta_full), value = TRUE))
            predicted_full = as.vector(X_full[, non_cell_terms, drop = FALSE] %*% beta_full[non_cell_terms])

            # Evaluating the interaction terms by coefficient name so rank-dropped columns are ignored safely.
            interaction_terms = grep(pattern = ":", x = names(beta_full), value = TRUE)
            interaction_terms = interaction_terms[grepl("group", interaction_terms)]
            if (length(interaction_terms) > 0) {
                interaction_df = X_full[, interaction_terms, drop = FALSE]
                interaction_coef = beta_full[interaction_terms]
                interaction_estimates_all = as.vector(interaction_df %*% interaction_coef)
            } else {
                interaction_estimates_all = rep(0, nrow(data))
            }
            group_terms = grep(pattern = "^group", x = names(beta_full), value = TRUE)
            if (length(group_terms) > 0) {
                interaction_estimates_all = interaction_estimates_all + beta_full[group_terms[1]]
            }

            # Likelihood Ratio Test
            vals <- lrtest(null, full)$Chisq[2]
            list(chisq = vals, predicted_full = predicted_full, interaction_estimates_all = interaction_estimates_all)
        }, error = function(e) {
            list(chisq = -1, predicted_full = rep(0, nrow(data)), interaction_estimates_all = rep(0, nrow(data)))
        })

        chisq = model_result$chisq
        generate_vector <- function(values, indices) {
            length_result <- length(values) + length(indices)
            result <- rep(0, length_result)
            value_positions <- setdiff(seq_along(result), indices)
            result[value_positions] <- values
            return(result)
        }    
        
        predictions = generate_vector(model_result$predicted_full, rm_index)
        interaction = generate_vector(model_result$interaction_estimates_all, rm_index)
    }
    """

    robjects.r(r_script)

    test_statistics = robjects.globalenv['chisq']
    test_statistics = np.array(test_statistics).astype(np.float64)

    predictions = robjects.globalenv['predictions']
    predictions = np.array(predictions).astype(np.float64)

    interaction = robjects.globalenv['interaction']
    interaction = np.array(interaction).astype(np.float64)
        
    return test_statistics, predictions, interaction

class calculate_test_statistic_multi_rep:
    def __init__(self, adata_list, group_id, k = 300, cell_type = None):
        self.adata_list = adata_list
        self.k = k
        self.group_id = group_id
        self.cell_type = cell_type
    
    def __call__(self, d, g):
        """
        Calculate the test statistic for a given gene and embedding dimension
        
        Inputs:
        g (int): The gene index
        d (int): The embedding dimension
        
        Returns:
        test_statistic (float): The likelihood ratio test statistic
        """
        
        # Extract the embedding and expression data
        emb = np.concatenate([self.adata_list[i].obsm['SpaceExpress'][:, d] for i in range(len(self.adata_list))])
        g_id = np.concatenate([np.array([self.group_id[i]]*self.adata_list[i].X.shape[0]) for i in range(len(self.adata_list))])
        rep = np.concatenate([np.array([i]*self.adata_list[i].X.shape[0]) for i in range(len(self.adata_list))])
        exp = np.concatenate([_as_dense_vector(self.adata_list[i].X[:, g]) for i in range(len(self.adata_list))])
        
        n_cells_total = sum(adata.X.shape[0] for adata in self.adata_list)
        if self.cell_type != None:
            ct = np.concatenate([self.adata_list[i].obs[self.cell_type].values.tolist() for i in range(len(self.adata_list))])
            df = pd.DataFrame({'emb': emb, 'exp': exp, 'group_id': g_id, 'rep': rep, 'cell_type':ct})
            
            try:    
                test_statistic = spline_multi_rep_ct(df, self.k)
                return test_statistic
            except Exception as e:
                warnings.warn(f"DSE fit failed at dim={d}, gene={g}: {e}", RuntimeWarning)
                return _fallback_dse_result(n_cells_total)
            
        else:
            df = pd.DataFrame({'emb': emb, 'exp': exp, 'group_id': g_id, 'rep': rep})
            
            try:    
                test_statistic = spline_multi_rep(df, self.k)
                return test_statistic
            except Exception as e:
                warnings.warn(f"DSE fit failed at dim={d}, gene={g}: {e}", RuntimeWarning)
                return _fallback_dse_result(n_cells_total)


class calculate_test_statistic:
    def __init__(self, adata_1, adata_2, k = 300):
        self.adata_1 = adata_1
        self.adata_2 = adata_2
        self.k = k
    
    def __call__(self, d, g):
        """
        Calculate the test statistic for a given gene and embedding dimension
        
        Inputs:
        d (int): The embedding dimension
        g (int): The gene index
        
        Returns:
        test_statistic (float): The likelihood ratio test statistic
        """
        n_cells_total = self.adata_1.X.shape[0] + self.adata_2.X.shape[0]
        # Extract the embedding and expression data
        emb = np.concatenate([self.adata_1.obsm['SpaceExpress'][:, d], self.adata_2.obsm['SpaceExpress'][:, d]])
        data_id = np.concatenate([np.array([0]*self.adata_1.X.shape[0]), np.array([1]*self.adata_2.X.shape[0])])
        exp = np.concatenate([_as_dense_vector(self.adata_1.X[:, g]), _as_dense_vector(self.adata_2.X[:, g])])
        df = pd.DataFrame({'emb': emb, 'exp': exp, 'data_id': data_id})
        # Calculate the test statistic
        try:    
            test_statistic = spline(df, self.k)
            return test_statistic
        except Exception as e:
            warnings.warn(f"DSE fit failed at dim={d}, gene={g}: {e}", RuntimeWarning)
            return _fallback_dse_result(n_cells_total)
        

def empirical_null(df, quant_val = 0.75):
    """
    Adjusting the test statistics for the empirical null distribution
    
    Parameters:
    df (pd.DataFrame): A pandas DataFrame containing the test statistics
    
    Returns:
    df_fdr (pd.DataFrame): A pandas DataFrame containing the false discovery rates
    """
    # Convert pandas DataFrame to R DataFrame
    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_df = robjects.conversion.py2rpy(df)

    # Assign the DataFrame to an R variable
    robjects.globalenv['df'] = r_df
    robjects.globalenv['quant_val'] = quant_val

    # Define the R script as a Python string
    r_script = """
    library(fitdistrplus)
    
    n = ncol(df)
    quant_val = quant_val

    fdr = matrix(NA, nrow = nrow(df), ncol = ncol(df))
    for (i in 1:nrow(df)){
        # Extract the test statistics of embedding dimension i
        T <- as.numeric(df[i,])
        rm_index <- which(T == -1)
        T <- T[!T %in% T[rm_index]]
        T <- T[is.finite(T)]

        generate_vector <- function(values, indices) {
            length_result <- length(values) + length(indices)
            result <- rep(0, length_result)
            value_positions <- setdiff(seq_along(result), indices)
            result[value_positions] <- values
            return(result)
        }

        if (length(T) < 2 || var(T) == 0) {
            fdr[i,] = generate_vector(rep(1, length(T)), rm_index)
            next
        }
        
        # Define the quantile value for thresholding
        q = quantile(T, quant_val)
        
        # Extract the subset of T values below the quantile threshold
        A0 = T[T < q]
        A0 = A0[is.finite(A0) & A0 > 0]
        if (length(A0) < 2 || var(A0) == 0) {
            fdr[i,] = generate_vector(rep(1, length(T)), rm_index)
            next
        }
        
        # Fit a gamma distribution to the subset of T values
        fit_a0 = tryCatch(
            fitdist(data = A0, distr = "gamma", method = "mle"),
            error = function(e) NULL
        )
        if (is.null(fit_a0)) {
            fdr[i,] = generate_vector(rep(1, length(T)), rm_index)
            next
        }
        
        # Extract the shape (k) and rate (1/theta) parameters of the fitted gamma distribution
        k = fit_a0$estimate[1]
        theta = 1 / fit_a0$estimate[2]
        
        # Calculate the proportion of T values below the quantile threshold
        n0 = sum(T < q)
        p0 = n0 / (n * pgamma(q, shape = k, scale = theta))
        
        # Calculate the false discovery rate (FDR) 
        num = (1 - pgamma(T, shape = k, scale = theta))
        denom = 1 - ecdf(T)(T)
        fdr_new = p0 * num / (denom + 1e-5)
        
        fdr[i,] = generate_vector(fdr_new, rm_index)
    } 
    """

    # Execute the R script
    robjects.r(r_script)

    fdr = robjects.globalenv['fdr']
    fdr = np.array(fdr).astype(np.float64)
    df_fdr = pd.DataFrame(fdr, columns = df.columns, index = df.index)
    
    return df_fdr


def SpaceExpress_DSE(emb, adata_list, cell_type = None, k = 300, n_jobs=-1, multi = False, 
                     group_id = None, quant_val = 0.75):
    """
    Perform the SpaceExpress differential spatial expression analysis
    
    Inputs:
    adata_list (list): A list of anndata objects containing the spatial expression data
    n_jobs (int): The number of parallel jobs to run
    
    Returns:
    df_fdr (pd.DataFrame): A pandas DataFrame containing the false discovery rates
    """
    _warn_if_k_is_large(adata_list, k, multi=multi, cell_type=cell_type, group_id=group_id)

    for i in range(len(adata_list)):
        adata_list[i].obsm['SpaceExpress'] = emb[i]
        
    # Extract the data from the anndata objects
    adata_1 = adata_list[0]
    adata_2 = adata_list[1]

    num_gene = adata_1.X.shape[1]
    num_dim = adata_1.obsm['SpaceExpress'].shape[1]    

    # Create a function to calculate the test statistic
    
    if multi == True:
        assert group_id != None, "There is no group_id"
        if len(group_id) != len(adata_list):
            raise ValueError("group_id length must match adata_list length.")
        cal_statistic = calculate_test_statistic_multi_rep(adata_list, cell_type = cell_type, k = k, group_id = group_id)
    else:    
        cal_statistic = calculate_test_statistic(adata_1, adata_2, k = k)

    # Create a list of jobs
    jobs = [(d, g) for d in range(num_dim) for g in range(num_gene)]

    # Run jobs in parallel with progress monitoring
    results = Parallel(n_jobs=n_jobs)(delayed(cal_statistic)(d, g) for d, g in tqdm(jobs, desc="Processing"))

    # Convert results back to the test_statistics array
    test_statistics = np.zeros((num_dim, num_gene))
    # predictions_1 = np.zeros((adata_1.shape[0], num_gene, num_dim))
    # predictions_2 = np.zeros((adata_2.shape[0], num_gene, num_dim))
    # interaction_1 = np.zeros((adata_1.shape[0], num_gene, num_dim))
    # interaction_2 = np.zeros((adata_2.shape[0], num_gene, num_dim))
    
    num_cell_list = [adata.shape[0] for adata in adata_list]
    predictions_list, interactions_list = [np.zeros((adata.shape[0], num_gene, num_dim)) for adata in adata_list], [np.zeros((adata.shape[0], num_gene, num_dim)) for adata in adata_list]
    for idx, (d, g) in enumerate(jobs):
        test_stat, pred, inter = results[idx]
        # print(f"test_stat: {test_stat}, shape: {test_stat[0].shape}, test_statistics: {test_statistics.shape}")
        test_statistics[d, g] = test_stat[0]
        
        idx1, idx2 = 0, 0
        for i in range(len(adata_list)):
            idx1 = idx2
            idx2 += num_cell_list[i]
            predictions_list[i][:, g, d] = pred[idx1:idx2]
            interactions_list[i][:, g, d] = inter[idx1:idx2]

            # predictions_list[i][:, g, d] = pred[num_cell_list[i]*i:num_cell_list[i]*(i+1)]
            # interactions_list[i][:, g, d] = inter[num_cell_list[i]*i:num_cell_list[i]*(i+1)]
            
    df_test_statistics = pd.DataFrame(data=test_statistics, index=range(num_dim), columns=list(adata_list[0].var_names))
    df_fdr = empirical_null(df_test_statistics, quant_val)

    for i in range(len(adata_list)):
        adata_list[i].varm['DSE-fdr'] = df_fdr.T
    
    for i in range(len(adata_list)):
        adata_list[i].obsm['DSE-pred'] = predictions_list[i]
        adata_list[i].obsm['DSE-inter'] = interactions_list[i]
    # adata_list[0].obsm['DSE-pred'] = predictions_1
    # adata_list[1].obsm['DSE-pred'] = predictions_2
    # adata_list[0].obsm['DSE-inter'] = interaction_1
    # adata_list[1].obsm['DSE-inter'] = interaction_2 

    return df_fdr, adata_list

def summary_DSE (df, threshold = 0.001):
    """
    Summarize the differential spatial expression results
    
    Inputs:
    df (pd.DataFrame): A pandas DataFrame containing the false discovery rates
    
    Returns:
    summary_df (pd.DataFrame): A pandas DataFrame containing the summary results
    """
    # Identify the positions where the condition is met
    condition = df < threshold

    # Create a dictionary to store the summary data
    summary_dict = {'Dim': [], 'DSE Count': [], 'DSE': []}

    # Iterate over each row to summarize the data
    for row_index, row_data in condition.iterrows():
        selected_columns = row_data[row_data].index.tolist()
        if selected_columns:
            summary_dict['Dim'].append(row_index)
            summary_dict['DSE Count'].append(len(selected_columns))
            summary_dict['DSE'].append(selected_columns)
        else:
            summary_dict['Dim'].append(row_index)
            summary_dict['DSE Count'].append(0)
            summary_dict['DSE'].append([])

    # Convert the summary dictionary to a DataFrame
    summary_df = pd.DataFrame(summary_dict)
    return summary_df
