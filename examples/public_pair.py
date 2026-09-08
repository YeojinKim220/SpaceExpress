#!/usr/bin/env python3
"""Reproducible two-sample run; paths in a config are relative to that config."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
import pickle
import subprocess
import sys
import time
import traceback
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import torch
from scipy.sparse.csgraph import connected_components
from sklearn.neighbors import kneighbors_graph

import SpaceExpress as se


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    temporary.replace(path)


def balanced_positions(frame, n, seed):
    if len(frame) < n:
        raise ValueError(f'Requested {n} observations but only {len(frame)} pass QC')
    if len(frame) == n:
        return frame.index.to_numpy()
    work = frame.copy()
    for axis in ('x', 'y'):
        values = work[axis].to_numpy(float)
        edges = np.linspace(values.min(), values.max(), 11)
        work[f'_{axis}bin'] = np.clip(np.digitize(values, edges[1:-1]), 0, 9)
    rng = np.random.default_rng(seed)
    selected = []
    for indices in work.groupby(['_xbin', '_ybin'], observed=True).groups.values():
        quota = int(np.floor(n * len(indices) / len(work)))
        selected.extend(rng.choice(np.asarray(indices), quota, replace=False).tolist())
    pool = work.index[~work.index.isin(selected)].to_numpy()
    selected.extend(rng.choice(pool, n-len(selected), replace=False).tolist())
    return np.asarray(selected)


def count_adata(matrix, obs, genes, coords):
    matrix = sp.csr_matrix(matrix, dtype=np.float32)
    if not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
        raise ValueError('Raw counts must be finite and nonnegative')
    if not np.allclose(matrix.data, np.rint(matrix.data)):
        raise ValueError('Expected raw integer counts, not normalized expression')
    result = ad.AnnData(matrix, obs=obs.copy(), var=pd.DataFrame(index=genes))
    if not result.obs_names.is_unique or not result.var_names.is_unique:
        raise ValueError('Duplicate observation or gene names in raw input')
    result.obsm['spatial'] = np.asarray(coords, dtype=np.float32)
    result.obs['total_counts'] = np.asarray(matrix.sum(axis=1)).ravel()
    result.obs['n_genes_by_counts'] = np.asarray((matrix > 0).sum(axis=1)).ravel()
    mt = result.var_names.str.upper().str.startswith('MT-')
    result.obs['pct_counts_mt'] = 100 * np.asarray(matrix[:, mt].sum(axis=1)).ravel() / np.maximum(result.obs['total_counts'], 1)
    return result


def load_slide(sample):
    locations = pd.read_csv(sample['locations']).set_index('barcode')
    labels = pd.read_csv(sample['cell_types']).drop(columns=['Unnamed: 0'], errors='ignore').set_index('barcode')
    if not locations.index.is_unique or not labels.index.is_unique:
        raise ValueError('Duplicate barcodes in locations or cell types')
    with open(sample['dge'], newline='') as handle:
        genes = next(csv.reader(handle))[1:]
    matrices, barcodes = [], []
    for chunk in pd.read_csv(sample['dge'], chunksize=256,
                             dtype={g: np.float32 for g in genes}):
        barcodes.extend(chunk['barcode'].astype(str))
        matrices.append(sp.csr_matrix(chunk[genes].to_numpy(dtype=np.float32)))
    matrix = sp.vstack(matrices, format='csr')
    obs = pd.DataFrame(index=pd.Index(barcodes, name='barcode'))
    obs = obs.join(locations[['x', 'y']]).join(labels[['max_cell_type']])
    obs['annotation'] = obs['max_cell_type'].astype('string').fillna('Unknown').astype(str)
    return count_adata(matrix, obs, genes, obs[['x', 'y']].to_numpy())


def load_stereo(sample, root, index):
    backed = ad.read_h5ad(sample['h5ad'], backed='r')
    try:
        full_obs = backed.obs.copy()
        full_obs[['x', 'y']] = np.asarray(backed.obsm['spatial'])[:, :2]
        full_obs.to_csv(root / f'source_observations_{index}.csv')
        mask = full_obs['annotation'].astype(str).eq(sample['annotation']).to_numpy()
        sub = backed[np.flatnonzero(mask), :].to_memory()
        source_count = backed.n_obs
    finally:
        backed.file.close()
    counts = sub.layers['count']
    result = count_adata(counts, sub.obs, sub.var_names, sub.obsm['spatial'])
    return result, source_count


def graph_policy(pair):
    scans, choices = [], []
    for a in pair:
        counts = {}
        for k in range(4, 21):
            graph = kneighbors_graph(a.obsm['spatial'], k, include_self=False)
            counts[k] = int(connected_components(graph.maximum(graph.T), directed=False)[0])
        choices.append(min(counts, key=lambda k: (counts[k], k)))
        scans.append(counts)
    k = max(choices)
    return k, [{'components_by_k': scan, 'components_at_selected_k': scan[k]} for scan in scans]


def prepare(config, root, status):
    pair, summaries = [], []
    for i, (sample, target) in enumerate(zip(config['samples'], config['target_spots'])):
        if config['type'] == 'slide':
            a = load_slide(sample)
            source_count = a.n_obs
        else:
            a, source_count = load_stereo(sample, root, i)
        selected_tissue_count = a.n_obs
        coords = a.obsm['spatial']
        eligible = np.isfinite(coords).all(axis=1) & (a.obs['n_genes_by_counts'].to_numpy() >= 3)
        frame = a.obs.copy()
        frame[['x', 'y']] = coords
        frame['qc_pass'] = eligible
        candidates = pd.DataFrame(coords[eligible], index=np.flatnonzero(eligible), columns=['x', 'y'])
        chosen = np.sort(balanced_positions(candidates, target, config['seed'] + i))
        frame['selected_for_training'] = False
        frame.iloc[chosen, frame.columns.get_loc('selected_for_training')] = True
        frame.to_csv(root / f'qc_observations_{i}.csv')
        a = a[chosen].copy()
        a.obs['source_sample'] = sample['name']
        a.obs['condition'] = sample['condition']
        a.layers['counts'] = a.X.copy()
        a.write_h5ad(root / f'raw_selected_{i}.h5ad', compression='gzip')
        summaries.append({'sample': sample['name'], 'original_observations': source_count,
                          'tissue_observations': selected_tissue_count,
                          'qc_observations': int(eligible.sum()), 'training_observations': a.n_obs,
                          'original_genes': a.n_vars})
        pair.append(a)
    common = sorted(set(pair[0].var_names).intersection(pair[1].var_names))
    pair = [a[:, common].copy() for a in pair]
    keep = np.logical_and.reduce([np.asarray((a.X > 0).sum(axis=0)).ravel() >= 3 for a in pair])
    pair = [a[:, keep].copy() for a in pair]
    for a in pair:
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    pair, hvg, diagnostics = se.select_hvg_after_outlier(pair, n_top_genes=config['n_hvg'], z_threshold=4)
    graph_k, graph_scan = graph_policy(pair)
    clipping = []
    for i, a in enumerate(pair):
        dense = a.X.toarray()
        q95 = np.percentile(dense, 95, axis=0)
        zeros = a.var_names[(q95 == 0) & (dense.max(axis=0) > 0)].tolist()
        clipping.append({'sample': config['samples'][i]['name'],
                         'genes_zeroed_by_embedding_q95': zeros,
                         'genes_already_all_zero_after_4sd': a.var_names[dense.max(axis=0) == 0].tolist()})
        pd.DataFrame({'gene': a.var_names, 'embedding_q95': q95,
                      'positive_observations_after_4sd': (dense > 0).sum(axis=0)}).to_csv(root / f'gene_qc_{i}.csv', index=False)
        if not np.isfinite(dense).all():
            raise ValueError('Nonfinite prepared expression')
        a.write_h5ad(root / f'prepared_condition{i}.h5ad', compression='gzip')
    pd.Series(hvg, name='gene').to_csv(root / 'hvg_genes.csv', index=False)
    pd.DataFrame(summaries).to_csv(root / 'sample_summary.csv', index=False)
    status.update(samples=summaries, hvg_preprocessing=diagnostics, clipping=clipping,
                  graph_k=graph_k, graph_scan=graph_scan)


def embed(config, root, status, probe=False):
    if not torch.cuda.is_available():
        raise RuntimeError('CUDA required for the full-size embedding stages')
    pair = [ad.read_h5ad(root / f'prepared_condition{i}.h5ad') for i in range(2)]
    if [a.n_obs for a in pair] != config['target_spots']:
        raise ValueError('Prepared observation counts do not match config')
    graph_k = json.loads((root / 'prepare_status.json').read_text())['graph_k']
    shortest = []
    for i, a in enumerate(pair):
        path = root / f'shortest_condition{i}.pkl'
        if not path.exists():
            temporary = path.with_suffix('.partial')
            if temporary.exists():
                raise RuntimeError(f'Incomplete shortest paths exist: {temporary}; inspect before retrying')
            se.shortest_path(a.obsm['spatial'], str(temporary), k=graph_k)
            temporary.replace(path)
        shortest.append(str(path))
    torch.cuda.reset_peak_memory_stats()
    start = time.perf_counter()
    epochs = 2 if probe else config['epochs']
    embeddings, model = se.train_SpaceExpress_multi(pair, shortest, epochs=epochs,
        patience=config['patience'], emb_dim=4, hid_dim=32, random_seed=config['seed'], save_model=True)
    status.update(device=torch.cuda.get_device_name(), epochs_requested=epochs,
                  training_seconds=time.perf_counter()-start,
                  peak_cuda_allocated_gib=torch.cuda.max_memory_allocated()/2**30,
                  peak_cuda_reserved_gib=torch.cuda.max_memory_reserved()/2**30,
                  device_memory_gib=torch.cuda.get_device_properties(0).total_memory/2**30)
    for a, emb in zip(pair, embeddings):
        if emb.shape != (a.n_obs, 4) or not np.isfinite(emb).all() or np.any(emb.std(axis=0) < 1e-8):
            raise ValueError('Nonfinite, collapsed, or incorrectly sized embedding')
    status['embedding_shapes'] = [list(e.shape) for e in embeddings]
    if not probe:
        for i, (a, emb) in enumerate(zip(pair, embeddings)):
            a.obsm['SpaceExpress'] = emb
            a.write_h5ad(root / f'embedded_condition{i}.h5ad', compression='gzip')
        with (root / 'embedding.pkl').open('wb') as handle:
            pickle.dump(embeddings, handle)
        torch.save(model.state_dict(), root / 'model_state.pt')


def dse(config, root, status, k):
    out = root / f'k{k}'
    out.mkdir(exist_ok=True)
    pair = [ad.read_h5ad(root / f'embedded_condition{i}.h5ad') for i in range(2)]
    embeddings = [a.obsm['SpaceExpress'] for a in pair]
    fdr, fitted = se.SpaceExpress_DSE(embeddings, pair, k=k,
        n_jobs=int(os.environ.get('SLURM_CPUS_PER_TASK', '4')), multi=False)
    values = fdr.to_numpy(float)
    failed = fitted[0].varm['DSE-fit-failed'].to_numpy().T
    if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
        raise ValueError('Invalid FDR range')
    if not np.all(values[failed] == 1):
        raise ValueError('Failed fits did not receive FDR=1')
    if any(a.uns['spaceexpress_dse']['remove_mean_sd_outliers'] for a in fitted):
        raise ValueError('Duplicate mean+4SD removal was unexpectedly enabled')
    if not all(np.isfinite(a.obsm[key]).all() for a in fitted for key in ['DSE-pred', 'DSE-inter']):
        raise ValueError('Nonfinite predictions')
    if failed.all():
        raise ValueError('All DSE fits failed')
    fdr.to_csv(out / 'fdr.csv')
    fitted[0].varm['DSE-statistic'].to_csv(out / 'test_statistics.csv')
    fitted[0].varm['DSE-fit-failed'].to_csv(out / 'fit_failed.csv')
    with (out / 'result.pkl').open('wb') as handle:
        pickle.dump([fdr, fitted], handle, protocol=pickle.HIGHEST_PROTOCOL)
    genes = fdr.columns[(fdr < 0.001).any(axis=0)].tolist()
    status.update(k=k, fdr_min=float(values.min()), fdr_max=float(values.max()),
                  failed_gene_dimensions=int(failed.sum()),
                  failed_genes=fdr.columns[failed.any(axis=0)].tolist(),
                  dse_gene_count_0_001=len(genes), dse_genes_0_001=genes,
                  repeated_mean_sd_removal=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--stage', choices=['prepare', 'probe', 'embed', 'dse', 'report'], required=True)
    parser.add_argument('--k', type=int, choices=[30, 50, 100], default=30)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = json.loads(config_path.read_text())
    for sample in config['samples']:
        for key in ('dge', 'locations', 'cell_types', 'h5ad'):
            if key in sample:
                sample[key] = str((config_path.parent / sample[key]).resolve())
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    stage = f'dse_k{args.k}' if args.stage == 'dse' else args.stage
    status_path = root / f'{stage}_status.json'
    if status_path.exists() and json.loads(status_path.read_text())['status'] == 'PASS':
        raise RuntimeError(f'{stage} already completed; use a new output directory for a fresh run')
    version = importlib.metadata.version('spaceexpress')
    if version != '0.1.5':
        raise RuntimeError(f'Expected installed release 0.1.5, found {version}')
    snapshot = root / 'config_resolved.json'
    if snapshot.exists():
        if json.loads(snapshot.read_text()) != config:
            raise ValueError('Config changed within a run; choose a new output directory')
    else:
        write_json(snapshot, config)
        (root / 'pip_freeze.txt').write_text(subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'], text=True))
    sources = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(se.__file__).parent.glob('*.py')}
    status = {'status': 'RUNNING', 'stage': stage, 'version': version,
              'python': sys.executable, 'package_path': se.__file__, 'source_sha256': sources,
              'slurm_job_id': os.environ.get('SLURM_JOB_ID')}
    write_json(status_path, status)
    start = time.perf_counter()
    try:
        if args.stage == 'prepare':
            prepare(config, root, status)
        elif args.stage in ('probe', 'embed'):
            embed(config, root, status, probe=args.stage == 'probe')
        elif args.stage == 'dse':
            dse(config, root, status, args.k)
        else:
            from public_pair_notebook import build_and_execute
            build_and_execute(root)
        status['status'] = 'PASS'
    except Exception as exc:
        status.update(status='ERROR', error=repr(exc), traceback=traceback.format_exc())
        raise
    finally:
        status['elapsed_seconds'] = time.perf_counter()-start
        write_json(status_path, status)


if __name__ == '__main__':
    main()
