"""Shared expression preprocessing before embedding and DSE."""

from __future__ import annotations

import warnings

import numpy as np
import scanpy as sc
import scipy.sparse as sp


def preprocessing(adata_list, n_top_genes=1000, z_threshold=4.0, flavor="seurat"):
    """Select common genes, zero pooled outliers, and take per-sample HVG union.

    Normalize and log-transform each sample before calling this function when
    using the default Seurat flavor. This function does not normalize, log,
    clip percentiles, or scale expression. For other input representations,
    callers must supply data appropriate for the chosen Scanpy HVG flavor.

    Values at or above the pooled gene mean + z_threshold * sample SD are
    set to zero, without deleting observations. Then each sample selects
    n_top_genes HVGs (default 1000) from the common genes, following Scanpy
    cutoff-tie behavior (which can retain more than the requested number). All genes selected
    in at least one sample are retained in every sample, in sorted order.
    If there are at most n_top_genes common genes, all are retained.

    Returns a list of independent AnnData copies with cleaned CSR expression;
    input objects are not modified. Other layers are only subsetted, so raw
    counts can be preserved separately. Per-sample HVG flags are in
    var['highly_variable']; provenance is in uns['spaceexpress_preprocessing'].
    Call once before passing the same list to training and DSE.
    """
    if len(adata_list) < 2:
        raise ValueError("At least two AnnData objects are required.")
    if isinstance(n_top_genes, (bool, np.bool_)) or not isinstance(n_top_genes, (int, np.integer)) or n_top_genes < 1:
        raise ValueError("n_top_genes must be a positive integer.")
    if not np.isfinite(z_threshold) or z_threshold <= 0:
        raise ValueError("z_threshold must be finite and positive.")
    for i, item in enumerate(adata_list):
        if not item.var_names.is_unique:
            raise ValueError(f"Sample {i} has duplicate gene names.")
        if item.n_obs < 2:
            raise ValueError(f"Sample {i} needs at least two observations.")
        values = item.X.data if sp.issparse(item.X) else np.asarray(item.X)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f"Sample {i} expression must be finite and nonnegative.")

    common = sorted(set.intersection(*(set(item.var_names) for item in adata_list)))
    if not common:
        raise ValueError("No common genes across samples.")
    cleaned = [item[:, common].copy() for item in adata_list]
    matrices = [sp.csr_matrix(item.X, dtype=np.float64, copy=True) for item in cleaned]
    for matrix in matrices:
        matrix.sum_duplicates()
        matrix.eliminate_zeros()

    pooled = sp.vstack(matrices, format="csr")
    n_obs = pooled.shape[0]
    means = np.asarray(pooled.sum(axis=0)).ravel() / n_obs
    sums_of_squares = np.asarray(pooled.power(2).sum(axis=0)).ravel()
    variances = (sums_of_squares - n_obs * means**2) / (n_obs - 1)
    thresholds = means + z_threshold * np.sqrt(np.maximum(variances, 0.0))

    removed_entries, hvg_sets = [], []
    for item, matrix in zip(cleaned, matrices):
        remove = matrix.data >= thresholds[matrix.indices]
        removed_entries.append(int(remove.sum()))
        matrix.data[remove] = 0.0
        matrix.eliminate_zeros()
        item.X = matrix.astype(np.float32)
        if len(common) > n_top_genes:
            sc.pp.highly_variable_genes(item, n_top_genes=n_top_genes, flavor=flavor)
        else:
            item.var["highly_variable"] = True
        hvg_sets.append(set(item.var_names[item.var["highly_variable"]]))

    hvg = sorted(set.union(*hvg_sets))
    if not hvg:
        raise ValueError("No HVGs selected from the cleaned expression.")
    result = [item[:, hvg].copy() for item in cleaned]
    for i, item in enumerate(result):
        item.uns["spaceexpress_preprocessing"] = {
            "mean_sd_outliers_removed": True,
            "method": "pooled_mean_plus_sample_sd",
            "z_threshold": float(z_threshold),
            "stage": "before_hvg",
            "hvg_selection": "per_sample_union",
            "hvg_flavor": flavor,
            "n_top_genes_per_sample": int(n_top_genes),
            "common_genes_before_hvg": len(common),
            "selected_hvg_per_sample": np.array([len(s) for s in hvg_sets]),
            "selected_hvg_union": len(hvg),
            "removed_expression_entries": np.array(removed_entries),
            "sample_index": i,
        }
    return result


def select_hvg_after_outlier(adata_list, n_top_genes=200, z_threshold=4.0):
    """Compatibility wrapper returning (samples, genes, diagnostics).

    Selection now uses the per-sample HVG union, so the result can contain more
    than n_top_genes genes. Prefer preprocessing(), whose default is 1000 per
    sample and whose return value is just the prepared list.
    """
    warnings.warn(
        "Use se.preprocessing(adata_list). select_hvg_after_outlier now uses "
        "per-sample HVG union rather than a fixed-size batch-aware selection.",
        DeprecationWarning, stacklevel=2,
    )
    result = preprocessing(adata_list, n_top_genes=n_top_genes, z_threshold=z_threshold)
    diagnostics = dict(result[0].uns["spaceexpress_preprocessing"])
    for key in ("selected_hvg_per_sample", "removed_expression_entries"):
        diagnostics[key] = diagnostics[key].tolist()
    diagnostics["selected_hvg_nonzero_spots"] = []
    for item in result:
        counts = np.asarray((item.X > 0).sum(axis=0)).ravel()
        diagnostics["selected_hvg_nonzero_spots"].append({
            "min": int(counts.min()), "q25": float(np.quantile(counts, .25)),
            "median": float(np.median(counts)), "q75": float(np.quantile(counts, .75)),
            "max": int(counts.max()),
        })
    return result, result[0].var_names.tolist(), diagnostics
