# SpaceExpress 0.1.5

Korean version: [RELEASE_0.1.5.md](RELEASE_0.1.5.md)

## Changes

- Failed gene-embedding-dimension fits (negative or nonfinite statistics) receive FDR=1. Dimensions without a valid empirical-null fit also return FDR=1.
- Final FDR values are bounded to [0, 1]. The existing empirical-null estimator is not replaced with another adjustment method.
- `select_hvg_after_outlier` operates on normalized, log-transformed common genes, setting values at or above the pooled mean + 4 sample standard deviations to 0 before batch-aware HVG selection. It modifies copies, not input objects.
- Initial preprocessing history is stored in `.uns['spaceexpress_preprocessing']`. DSE skips repeated mean+4SD removal when every input has this history. Unmarked legacy inputs retain the previous behavior; mixed histories require an explicit setting.
- `SpaceExpress_DSE(..., remove_mean_sd_outliers=False)` also explicitly disables subsequent mean+4SD removal. Existing training-time 95% clipping is unchanged. The separate 99% filter in multi-replicate DSE is also retained.
- Returned AnnData `.varm` includes `DSE-statistic` and `DSE-fit-failed`, so failures can be inspected without inferring them from 0 predictions.

## Installation

Python 3.11 or newer and R are required. Install R packages `lmtest`, `fitdistrplus`, `dplyr`, and `lme4` first, and configure Python to find the R shared library. Python dependencies are declared in package metadata. Use the PyTorch installation appropriate for your CUDA environment.

```bash
Rscript -e 'install.packages(c("lmtest", "fitdistrplus", "dplyr", "lme4"), repos="https://cloud.r-project.org")'
python -m pip install 'spaceexpress[notebooks]==0.1.5'
```

`examples/install_release.sh BASE_PYTHON NEW_ENV` installs the PyPI release into a new venv sharing an existing scientific/R environment. Package download provenance is recorded in `pip_install_report.json`, with dependencies in `pip_freeze.txt`. This is not a fully isolated installation of all dependencies from scratch.

## Running Data

`examples/public_pair.py` resolves paths relative to its config and does not overwrite originals. It checks raw counts for at least 3 detected genes and finite spatial coordinates, then applies spatial-bin proportional sampling if needed. Common genes detected in at least 3 observations in each sample are normalized to 10,000, log1p-transformed, and processed as above to select 200 HVGs. Selected raw data and QC metrics are also saved.

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
