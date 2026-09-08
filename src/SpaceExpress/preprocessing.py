"""Preprocessing helpers for sparse spatial-expression inputs."""

from __future__ import annotations

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp


def select_hvg_after_outlier(
    adata_list,
    n_top_genes=200,
    z_threshold=4.0,
):
    """Remove pooled high-expression outliers before batch-aware HVG selection.

    The input AnnData objects must already be normalized and log-transformed.
    For every common gene, values greater than or equal to pooled
    ``mean + z_threshold * sample_sd`` are set to zero. Seurat HVGs are then
    selected from the cleaned matrices with each AnnData treated as a batch.

    Returns
    -------
    cleaned_hvg : list[AnnData]
        Copies restricted to the same sorted HVGs and containing cleaned values.
    hvg : list[str]
        Sorted selected gene names.
    diagnostics : dict
        Counts and nonzero-prevalence summaries for reproducibility.
    """
    if len(adata_list) < 2:
        raise ValueError("At least two AnnData objects are required.")
    if n_top_genes < 1:
        raise ValueError("n_top_genes must be positive.")
    if z_threshold <= 0:
        raise ValueError("z_threshold must be positive.")

    common = sorted(set.intersection(*(set(item.var_names) for item in adata_list)))
    if len(common) < n_top_genes:
        raise ValueError(
            f"Need at least {n_top_genes} common genes, found {len(common)}."
        )
    cleaned = [item[:, common].copy() for item in adata_list]
    matrices = [
        item.X.tocsr().astype(np.float64, copy=True)
        if sp.issparse(item.X)
        else sp.csr_matrix(np.asarray(item.X, dtype=np.float64))
        for item in cleaned
    ]

    pooled = sp.vstack(matrices, format="csr")
    n_obs = pooled.shape[0]
    means = np.asarray(pooled.sum(axis=0)).ravel() / n_obs
    sums_of_squares = np.asarray(pooled.power(2).sum(axis=0)).ravel()
    variances = (sums_of_squares - n_obs * means**2) / max(n_obs - 1, 1)
    thresholds = means + z_threshold * np.sqrt(np.maximum(variances, 0.0))

    removed_entries = []
    for item, matrix in zip(cleaned, matrices):
        remove = matrix.data >= thresholds[matrix.indices]
        removed_entries.append(int(remove.sum()))
        matrix.data[remove] = 0.0
        matrix.eliminate_zeros()
        item.X = matrix.astype(np.float32)

    combined = ad.concat(
        cleaned,
        label="_hvg_batch",
        keys=[f"sample{i}" for i in range(len(cleaned))],
        index_unique="-",
    )
    sc.pp.highly_variable_genes(
        combined,
        n_top_genes=n_top_genes,
        flavor="seurat",
        batch_key="_hvg_batch",
    )
    hvg = sorted(combined.var_names[combined.var["highly_variable"]].tolist())
    if len(hvg) != n_top_genes:
        raise RuntimeError(f"Expected {n_top_genes} HVGs, found {len(hvg)}.")

    result = [item[:, hvg].copy() for item in cleaned]
    prevalence = []
    for item in result:
        counts = np.asarray((item.X > 0).sum(axis=0)).ravel()
        prevalence.append(
            {
                "min": int(counts.min()),
                "q25": float(np.quantile(counts, 0.25)),
                "median": float(np.median(counts)),
                "q75": float(np.quantile(counts, 0.75)),
                "max": int(counts.max()),
            }
        )
    diagnostics = {
        "method": "pooled_mean_plus_sample_sd",
        "z_threshold": float(z_threshold),
        "common_genes_before_hvg": len(common),
        "removed_expression_entries": removed_entries,
        "selected_hvg_nonzero_spots": prevalence,
    }
    for item in result:
        item.uns["spaceexpress_preprocessing"] = {
            "mean_sd_outliers_removed": True,
            "z_threshold": float(z_threshold),
            "stage": "before_hvg",
        }
    return result, hvg, diagnostics
