# SpaceExpress 0.1.5

Korean version: [RELEASE_0.1.5.md](RELEASE_0.1.5.md)

## Changes

- Failed gene-embedding-dimension fits (negative or nonfinite statistics) receive FDR=1. Dimensions without a valid empirical-null fit also return FDR=1.
- Final FDR values are bounded to [0, 1]. The existing empirical-null estimator is not replaced with another adjustment method.
- `se.preprocessing(adata_list)` selects common genes, sets expression at or above the pooled mean + 4 sample standard deviations to 0, selects 1000 HVGs per sample, and takes their union. Users normalize and log-transform outside the function. It returns copies of the inputs.
- Preprocessing history, per-sample HVG counts, and union size are stored in `.uns['spaceexpress_preprocessing']`. Training does not reselect HVGs for these union inputs.
- The DSE `remove_mean_sd_outliers` argument and mean+4SD removal have been deleted. DSE does not repeat this step regardless of preprocessing history. Training-time 95% clipping and the multi-replicate DSE 99% filter remain.
- Returned AnnData `.varm` includes `DSE-statistic` and `DSE-fit-failed`, so failures can be inspected without inferring them from 0 predictions.

## Shared Preprocessing API and Duplication Review

These API changes are unreleased changes in the current working tree. The PyPI 0.1.5 installation command below alone does not install the new API.

```python
for adata in adata_list:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

adata_list = se.preprocessing(adata_list)
emb = se.train_SpaceExpress_multi(adata_list, shortest_file_path_list)
fdr, adata_list = se.SpaceExpress_DSE(emb, adata_list)
```

Normalization can be adapted to the data; the default `flavor="seurat"` expects log-transformed, nonnegative expression. You can change `n_top_genes=1000`, `z_threshold=4.0`, and `flavor="seurat"`. When there are at most 1000 common genes, all are retained. Otherwise, per-sample Scanpy HVG results are combined. Cutoff ties follow Scanpy selection behavior and can make the per-sample count exceed the requested count. Union size is not fixed.

The return value is one AnnData list. All samples share the same sorted gene order, and `.X` contains cleaned CSR matrices. Observations are not deleted; other layers are only subsetted by gene, allowing raw counts to be preserved separately. `.var['highly_variable']` records selection in that sample. A gene marked False in one sample remains if another sample marks it True.

| Stage | Current processing | Duplication review |
| --- | --- | --- |
| Before the call | User-defined normalization and log transformation | Not performed inside the new function |
| `preprocessing` | Common genes, pooled mean+4SD zeroing, per-sample HVGs and union | Call once before training and DSE |
| `train_SpaceExpress_multi` | Per-sample 95% clipping, mean centering and SD scaling, gene-order alignment | Skips HVG reselection for new prepared inputs; preserves union |
| Training without the new preprocessing | Existing per-sample HVG selection followed by intersection | Legacy path retained; results can differ from the new path |
| Single-sample `train_SpaceExpress` | 95% clipping and scaling | Skips HVG reselection for new prepared inputs |
| Two-sample DSE | Per-group SD scaling without centering | Mean+4SD removal deleted |
| Multi-replicate DSE | Excludes observations at or above the pooled 99% quantile from fitting, per-replicate SD scaling | Mean+4SD removal deleted; 99% filter retained for subsequent review |

Training clips and scales copies. Pass the same prepared list, with its pre-training `.X` preserved, to DSE. Mixing new union inputs with unprepared inputs raises a training error. The training `num_hvg` argument does not trigger reselection for new union inputs.

The old `select_hvg_after_outlier` remains as a compatibility wrapper returning `(samples, genes, diagnostics)`. It now also uses the per-sample union, so results differ from the former fixed-size selection. The wrapper retains its old default of 200; the new API defaults to 1000 per sample. Remove the deleted `remove_mean_sd_outliers` argument from calling code as well.

## Installation

Python 3.11 or newer and R are required. Install R packages `lmtest`, `fitdistrplus`, `dplyr`, and `lme4` first, and configure Python to find the R shared library. Python dependencies are declared in package metadata. Use the PyTorch installation appropriate for your CUDA environment.

```bash
Rscript -e 'install.packages(c("lmtest", "fitdistrplus", "dplyr", "lme4"), repos="https://cloud.r-project.org")'
python -m pip install 'spaceexpress[notebooks]==0.1.5'
```

`examples/install_release.sh BASE_PYTHON NEW_ENV` installs the PyPI release into a new venv sharing an existing scientific/R environment. Package download provenance is recorded in `pip_install_report.json`, with dependencies in `pip_freeze.txt`. This is not a fully isolated installation of all dependencies from scratch.

## Running Data

`examples/public_pair.py` resolves paths relative to its config and does not overwrite originals. It checks raw counts for at least 3 detected genes and finite spatial coordinates, then applies spatial-bin proportional sampling if needed. Common genes detected in at least 3 observations in each sample are normalized to 10,000, log1p-transformed, and processed by the new per-sample HVG union. The config's `n_hvg` is the per-sample selection count and defaults to 1000 when omitted. Selected raw data and QC metrics are also saved.

```bash
bash examples/submit_public_pair.sh /path/to/env/bin/python /path/to/config_v015.json /path/to/results_v015
```

Default Slurm requests are 12 CPU/128 GB for preparation, 1 H100 80 GB GPU with 8 CPU/256 GB RAM for training, 12 CPU/128 GB per DSE job, and 4 CPU/64 GB for the notebook. Each job has a 4-hour limit. Set `SE_ACCOUNT`, `SE_QOS`, `SE_CPU_PARTITION`, `SE_GPU_PARTITION`, and `SE_GPU_CONSTRAINT` to configure your cluster. The default CPU partition is `cpu-medium`. For A100, set `SE_GPU_PARTITION=gpu-a100 SE_GPU_CONSTRAINT=A100-80GB`. These are requests, not measured requirements.

Stages run in order: preparation, a 2-epoch GPU probe using all selected observations, main training, k=30/50/100 DSE, and a results notebook. Dependencies require successful prior stages. Main training uses `epochs` and `patience` from the config. Reusing a run directory is rejected; choose a new output path for a new experiment.

Outputs include `sample_summary.csv`, QC and gene diagnostics CSVs, selected-raw/prepared/embedded H5ADs, model state, FDR/statistic/failure-mask CSVs, DSE pickles, per-stage JSON/time/memory logs, and an executed `results.ipynb` with figures. The notebook shows original versus actual observation counts, spatial coverage, counts and detected genes, embeddings, failed fits, k sensitivity, and gene expression and fitted results. Only open pickles from trusted runs.

## Interpretation and Limits

95% clipping can erase a gene's entire training signal in a sample when a sparse gene's upper bound is 0. This can occur after initial 4SD removal and is recorded in `gene_qc_*.csv` and the notebook. This release retains both agreed preprocessing steps and does not automatically substitute another method.

Slide-seq observations may be beads and Stereo-seq observations may be spatial bins; observation counts are not necessarily individual-cell counts. Inputs with only Brain tissue annotations are not assigned inferred detailed cell types for visualization. No additional mitochondrial-percentage or detection thresholds are applied automatically; inspect distributions and sample-specific QC before deciding.

FDR-range, failure-handling, and observation-preservation checks establish numerical behavior, not biological validity. Review k stability and spatial-neighbor recall as well. A comparison with 1 biological sample per condition does not replace replicated population-level inference.
