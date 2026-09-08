"""Build and execute a results notebook with embedded figures and CSV diagnostics."""
import sys
from pathlib import Path

import nbformat
from nbclient import NotebookClient


def build_and_execute(root):
    root = Path(root).resolve()
    nb = nbformat.v4.new_notebook()
    code = nbformat.v4.new_code_cell
    md = nbformat.v4.new_markdown_cell
    nb.cells = [
        md('# SpaceExpress 0.1.5: pair results\n\n'
           'Initial pooled mean + 4 sample SD zeroing before HVG selection and embedding 95th-percentile clipping are retained. '
           'Repeated mean + 4 SD removal in DSE is disabled. Numerical checks do not establish biological validity. '
           'Observations may be beads or spatial bins, not necessarily individual cells.'),
        code("""%matplotlib inline
from pathlib import Path
import json, pickle
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display
from sklearn.neighbors import NearestNeighbors
ROOT = Path.cwd()
K = 30
config = json.loads((ROOT / 'config_resolved.json').read_text())
figures = ROOT / 'figures'
figures.mkdir(exist_ok=True)
plt.rcParams.update({'figure.dpi': 110, 'font.size': 10})
def finish(fig, name):
    fig.tight_layout()
    fig.savefig(figures / f'{name}.png', dpi=160, bbox_inches='tight')
    plt.show()
    plt.close(fig)
def spatial(ax, xy, values, title, **kwargs):
    artist = ax.scatter(xy[:, 0], xy[:, 1], c=values, s=2, rasterized=True, **kwargs)
    ax.set(title=title, xlabel='Spatial x', ylabel='Spatial y')
    ax.set_aspect('equal')
    return artist
display(pd.read_csv(ROOT / 'sample_summary.csv'))
display(pd.DataFrame([json.loads((ROOT / f'{stage}_status.json').read_text()) for stage in ['prepare', 'probe', 'embed']])[
    ['stage', 'status', 'version', 'elapsed_seconds', 'peak_cuda_allocated_gib']])
pair = [ad.read_h5ad(ROOT / f'embedded_condition{i}.h5ad') for i in range(2)]
assert [a.n_obs for a in pair] == config['target_spots']
assert all(a.obsm['SpaceExpress'].shape == (a.n_obs, 4) for a in pair)
qc = [pd.read_csv(ROOT / f'qc_observations_{i}.csv', index_col=0) for i in range(2)]
"""),
        md('## Raw Spatial Coverage and QC\n\nGrey points show available tissue observations; colored points show the exact observations retained for training. Brain annotation is a tissue label, not a brain cell-type annotation.'),
        code("""categories = sorted(set.union(*(set(frame['annotation'].astype(str)) for frame in qc)))
colors = dict(zip(categories, plt.get_cmap('tab20')(np.linspace(0, 1, max(len(categories), 1)))))
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
composition = []
for i, (a, frame, sample) in enumerate(zip(pair, qc, config['samples'])):
    chosen = frame['selected_for_training'].astype(bool)
    xy = frame[['x', 'y']].to_numpy()
    axes[i, 0].scatter(xy[:, 0], xy[:, 1], c='#d6d6d6', s=2, rasterized=True)
    for label in categories:
        mask = chosen & frame['annotation'].astype(str).eq(label)
        axes[i, 0].scatter(xy[mask, 0], xy[mask, 1], s=2, color=colors[label], label=label, rasterized=True)
    axes[i, 0].set_title(f"{sample['name']}: {chosen.sum():,} selected / {len(frame):,} tissue")
    axes[i, 0].set_aspect('equal')
    axes[i, 0].legend(markerscale=4, fontsize=7)
    for j, field in enumerate(['total_counts', 'n_genes_by_counts'], start=1):
        artist = spatial(axes[i, j], xy, frame[field], f"{sample['name']}: {field}", cmap='viridis')
        fig.colorbar(artist, ax=axes[i, j], shrink=0.6)
    counts = pd.crosstab(frame['annotation'], frame['selected_for_training'])
    counts['sample'] = sample['name']
    composition.append(counts)
finish(fig, 'raw_spatial_qc')
display(pd.concat(composition))
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
for ax, field in zip(axes, ['total_counts', 'n_genes_by_counts', 'pct_counts_mt']):
    for frame, sample in zip(qc, config['samples']):
        ax.hist(frame[field], bins=70, alpha=0.5, label=sample['name'])
    ax.set(xlabel=field, ylabel='Observations')
    ax.legend()
finish(fig, 'qc_distributions')
"""),
        md('## Preprocessing Diagnostics\n\nThe 95% cap can become zero for sufficiently sparse genes. Such genes contribute no varying expression to embedding training in that sample. They are listed explicitly; the agreed preprocessing is not changed automatically.'),
        code("""prep = json.loads((ROOT / 'prepare_status.json').read_text())
display(prep['hvg_preprocessing'])
display(pd.DataFrame(prep['clipping']))
display(pd.DataFrame(prep['graph_scan']))
for i, sample in enumerate(config['samples']):
    gene_qc = pd.read_csv(ROOT / f'gene_qc_{i}.csv')
    print(sample['name'], 'genes with embedding q95 = 0:', int((gene_qc.embedding_q95 == 0).sum()))
    display(gene_qc.sort_values('positive_observations_after_4sd').head(15))
"""),
        md('## Embedding Coverage and Spatial Structure'),
        code("""fig, axes = plt.subplots(2, 4, figsize=(20, 10))
combined = np.concatenate([a.obsm['SpaceExpress'] for a in pair])
for i, (a, sample) in enumerate(zip(pair, config['samples'])):
    for d in range(4):
        artist = spatial(axes[i, d], a.obsm['spatial'], a.obsm['SpaceExpress'][:, d],
                         f"{sample['name']}: dimension {d+1}", cmap='coolwarm',
                         vmin=combined[:, d].min(), vmax=combined[:, d].max())
        fig.colorbar(artist, ax=axes[i, d], shrink=0.6)
finish(fig, 'embedding_spatial_dimensions')
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for a, sample in zip(pair, config['samples']):
    e = a.obsm['SpaceExpress']
    for ax, (d1, d2) in zip(axes, [(0, 1), (2, 3)]):
        ax.scatter(e[:, d1], e[:, d2], s=2, alpha=0.35, label=sample['name'], rasterized=True)
        ax.set(xlabel=f'Dimension {d1+1}', ylabel=f'Dimension {d2+1}')
        ax.legend(markerscale=4)
finish(fig, 'embedding_sample_comparison')
neighbor_diagnostics = []
for a, sample in zip(pair, config['samples']):
    spatial_nn = NearestNeighbors(n_neighbors=10).fit(a.obsm['spatial']).kneighbors(return_distance=False)
    embedding_nn = NearestNeighbors(n_neighbors=10).fit(a.obsm['SpaceExpress']).kneighbors(return_distance=False)
    overlap = np.array([len(set(x).intersection(y))/10 for x, y in zip(spatial_nn, embedding_nn)])
    neighbor_diagnostics.append({'sample': sample['name'], 'observations': a.n_obs,
                                 'mean_spatial_10nn_recall': overlap.mean(),
                                 'random_reference': 10/(a.n_obs-1),
                                 'embedding_dimension_sd': a.obsm['SpaceExpress'].std(axis=0).tolist()})
display(pd.DataFrame(neighbor_diagnostics))
pd.DataFrame(neighbor_diagnostics).to_csv(ROOT / 'embedding_neighbor_diagnostics.csv', index=False)
"""),
        md('## FDR, Failed Fits, and k Sensitivity\n\nFailed gene-dimension fits receive FDR=1. DSE sets use any dimension below 0.001. k sensitivity is descriptive; each condition has only one biological sample, so this is not a replicated population-level comparison.'),
        code("""statuses, sets = [], {}
for k in [30, 50, 100]:
    status = json.loads((ROOT / f'dse_k{k}_status.json').read_text())
    assert status['status'] == 'PASS'
    table = pd.read_csv(ROOT / f'k{k}' / 'fdr.csv', index_col=0)
    failed = pd.read_csv(ROOT / f'k{k}' / 'fit_failed.csv', index_col=0).to_numpy(bool).T
    assert np.isfinite(table.to_numpy()).all()
    assert table.to_numpy().min() >= 0 and table.to_numpy().max() <= 1
    assert np.all(table.to_numpy()[failed] == 1)
    assert not status['repeated_mean_sd_removal']
    statuses.append(status)
    sets[k] = set(table.columns[(table < 0.001).any(axis=0)])
display(pd.DataFrame(statuses)[['k', 'status', 'elapsed_seconds', 'fdr_min', 'fdr_max', 'failed_gene_dimensions', 'dse_gene_count_0_001']])
display(pd.DataFrame({f'k{k}': pd.Series(sorted(genes), dtype=str) for k, genes in sets.items()}))
print('Genes shared across all k:', sorted(set.intersection(*sets.values())))
with (ROOT / f'k{K}' / 'result.pkl').open('rb') as handle:
    fdr, fitted = pickle.load(handle)
ranking = pd.DataFrame({'min_fdr': fdr.min(axis=0),
                        'failed_dimensions': fitted[0].varm['DSE-fit-failed'].sum(axis=1)})
ranking = ranking.sort_values('min_fdr')
display(ranking.head(30))
display(ranking.loc[ranking.failed_dimensions > 0])
ranking.to_csv(ROOT / f'gene_ranking_k{K}.csv')
fig, ax = plt.subplots(figsize=(7, 4))
ax.hist(fdr.to_numpy().ravel(), bins=50)
ax.set(xlabel='FDR', ylabel='Gene-dimension tests', title=f'k = {K}')
finish(fig, f'fdr_distribution_k{K}')
"""),
        md('## Gene Expression and Fitted Spatial Patterns\n\nTop-ranked genes are shown even when none pass the FDR threshold. Expression is normalized/log-transformed and retains the initial 4SD zeroing. Prediction and interaction values are the package outputs.'),
        code("""for gene in ranking.index[:6]:
    g = fitted[0].var_names.get_loc(gene)
    d = int(fdr[gene].idxmin())
    fields = [[np.asarray(a.X[:, g].toarray()).ravel() for a in fitted],
              [a.obsm['DSE-pred'][:, g, d] for a in fitted],
              [a.obsm['DSE-inter'][:, g, d] for a in fitted]]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for j, (title, values) in enumerate(zip(['Observed expression', 'Fitted expression', 'Interaction'], fields)):
        lo, hi = min(v.min() for v in values), max(v.max() for v in values)
        for i, (a, sample) in enumerate(zip(fitted, config['samples'])):
            artist = spatial(axes[i, j], a.obsm['spatial'], values[i],
                             f"{sample['name']}: {title}", cmap='viridis', vmin=lo, vmax=hi)
            fig.colorbar(artist, ax=axes[i, j], shrink=0.6)
    fig.suptitle(f'{gene}: k={K}, dimension={d+1}, FDR={fdr.loc[d, gene]:.3g}', y=1.02)
    finish(fig, f'gene_{gene}_k{K}')
"""),
    ]
    nb.metadata['kernelspec'] = {'display_name': 'Python 3 (SpaceExpress 0.1.5)',
                                'language': 'python', 'name': 'python3'}
    path = root / 'results.ipynb'
    nbformat.write(nb, path)
    client = NotebookClient(nb, timeout=1800, kernel_name='python3',
                            resources={'metadata': {'path': str(root)}})
    client.km = client.create_kernel_manager()
    client.km.kernel_spec.argv = [sys.executable, '-m', 'ipykernel_launcher', '-f', '{connection_file}']
    try:
        client.execute()
    finally:
        nbformat.write(nb, path)
