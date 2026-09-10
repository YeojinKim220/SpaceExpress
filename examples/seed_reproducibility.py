"""Small independent-process GPU reproducibility check and executed notebook."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'src'
DATASETS = {
    'slide': 'SlideSeq_WT1_vs_Diabetes1',
    'brain': 'StereoSeq_E15.5_vs_E16.5_Brain',
}


def worker(args):
    if args.implementation == 'local':
        sys.path.insert(0, str(SOURCE))
    import anndata as ad
    import numpy as np
    import torch
    import SpaceExpress as se

    package = Path(se.__file__).resolve()
    if (SOURCE in package.parents) != (args.implementation == 'local'):
        raise RuntimeError(f'Unexpected package implementation: {package}')
    if args.device == 'cuda' and not torch.cuda.is_available():
        raise RuntimeError('GPU requested but unavailable; no silent CPU fallback')
    pair = [ad.read_h5ad(args.input / f'condition{i}.h5ad') for i in range(2)]
    kwargs = dict(device=args.device, epochs=args.epochs, patience=args.epochs + 1,
                  random_seed=args.seed, batch_size=32, save_model=True)
    if args.implementation == 'local':
        kwargs['deterministic'] = True
    started = time.perf_counter()
    if args.device == 'cuda':
        torch.cuda.reset_peak_memory_stats()
    embeddings, model = se.train_SpaceExpress_multi(pair,
        [str(args.input / f'shortest{i}.pkl') for i in range(2)], **kwargs)
    arrays = {f'embedding{i}': emb for i, emb in enumerate(embeddings)}
    arrays.update({f'weight_{key}': tensor.detach().cpu().numpy()
                   for key, tensor in model.state_dict().items()})
    if args.dse:
        prevalence = sum(np.asarray((a.X > 0).sum(axis=0)).ravel() for a in pair)
        genes = pair[0].var_names[np.argsort(-prevalence, kind='stable')[:8]]
        fdr, fitted = se.SpaceExpress_DSE(embeddings,
            [a[:, genes].copy() for a in pair], k=5, n_jobs=1, multi=False)
        arrays['fdr'] = fdr.to_numpy()
        arrays['dse_statistics'] = fitted[0].varm['DSE-statistic'].to_numpy()
        arrays['dse_failed'] = fitted[0].varm['DSE-fit-failed'].to_numpy()
        for i, a in enumerate(fitted):
            arrays[f'dse_predictions{i}'] = a.obsm['DSE-pred']
    assert all(np.isfinite(values).all() for values in arrays.values())
    assert all(np.all(emb.std(axis=0) > 1e-8) for emb in embeddings)
    np.savez(args.output.with_suffix('.npz'), **arrays)
    metadata = dict(seed=args.seed, hash_seed=os.environ.get('PYTHONHASHSEED'),
        implementation=args.implementation, package_path=str(package), torch=torch.__version__,
        source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in package.parent.glob('*.py')},
        device=torch.cuda.get_device_name() if args.device == 'cuda' else 'cpu',
        deterministic=torch.are_deterministic_algorithms_enabled(),
        epochs=args.epochs, shapes=[list(emb.shape) for emb in embeddings],
        seconds=time.perf_counter()-started,
        peak_cuda_allocated_gib=torch.cuda.max_memory_allocated()/2**30 if args.device == 'cuda' else 0)
    args.output.with_suffix('.json').write_text(json.dumps(metadata, indent=2) + '\n')


def compare(root, dataset, left, right, expected):
    import numpy as np
    with np.load(root / dataset / f'{left}.npz') as a, np.load(root / dataset / f'{right}.npz') as b:
        keys = [key for key in a.files if key.startswith(('embedding', 'weight_'))]
        exact = all(np.array_equal(a[key], b[key]) for key in keys)
        delta = max(float(np.abs(a[f'embedding{i}'] - b[f'embedding{i}']).max()) for i in range(2))
        dse_keys = [key for key in a.files if key.startswith(('fdr', 'dse_'))]
        dse_equal = all(np.array_equal(a[key], b[key]) for key in dse_keys) if dse_keys and all(key in b.files for key in dse_keys) else None
    return dict(dataset=dataset, comparison=f'{left} vs {right}', exact_embedding_and_weights=exact,
                max_embedding_abs_difference=delta, exact_dse=dse_equal,
                expected=expected, passed=(exact == (expected == 'same')) if expected != 'observe' else True)


def notebook(root):
    import nbformat
    from nbclient import NotebookClient
    nb = nbformat.v4.new_notebook()
    nb.cells = [
        nbformat.v4.new_markdown_cell('# SpaceExpress: seed reproducibility smoke test\n\n'
            '512 observations per sample, 200 previously selected HVGs, 30 epochs by default. '
            'Each run starts a fresh Python process on the same GPU. Fixed A/B also compare DSE '
            '(8 genes, k=5); this checks repeatability, not biological accuracy or full-size convergence. '
            'The installed 0.1.5 baseline is kept unchanged; fixed runs explicitly import local source. '
            'PYTHONHASHSEED must be set before Python starts. Different devices/software versions may differ.'),
        nbformat.v4.new_code_cell("""%matplotlib inline
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display
ROOT = Path.cwd()
summary = json.loads((ROOT / 'summary.json').read_text())
print('Overall:', summary['status'], '| total seconds:', round(summary['elapsed_seconds'], 1))
display(pd.DataFrame(summary['comparisons']))
display(pd.DataFrame(summary['runs'])[['dataset', 'run', 'device', 'seed', 'hash_seed', 'epochs', 'seconds', 'peak_cuda_allocated_gib', 'package_path']])
assert summary['status'] == 'PASS'
"""),
        nbformat.v4.new_markdown_cell('## Independent-process comparisons\n\n'
            'Fixed A/B: same RNG and hash seeds. Fixed A/hash-changed: only Python hash seed changes. '
            'Fixed A/seed-changed: only RNG seed changes. Legacy comparison illustrates the old hash-order sensitivity.'),
        nbformat.v4.new_code_cell("""for dataset in ['slide', 'brain']:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    with np.load(ROOT / dataset / 'fixed_a.npz') as base:
        x = np.concatenate([base[f'embedding{i}'] for i in range(2)]).ravel()
    for ax, label in zip(axes, ['fixed_b', 'hash_changed', 'seed_changed']):
        with np.load(ROOT / dataset / f'{label}.npz') as other:
            y = np.concatenate([other[f'embedding{i}'] for i in range(2)]).ravel()
        ax.scatter(x, y, s=2, alpha=0.3, rasterized=True)
        lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
        ax.plot([lo, hi], [lo, hi], color='black', linewidth=1)
        ax.set(title=f'{dataset}: {label}', xlabel='Fixed A embedding values', ylabel='Comparison embedding values')
    fig.tight_layout()
    fig.savefig(ROOT / f'{dataset}_seed_comparison.png', dpi=150)
    plt.show()
    plt.close(fig)
"""),
        nbformat.v4.new_markdown_cell('## Running again\n\n'
            'Use a new output directory with `run_seed_reproducibility.sh --output /path/to/new_run`. '
            'The launcher sets PYTHONHASHSEED, CUBLAS_WORKSPACE_CONFIG and thread limits before starting Python, '
            'and uses nohup with a 25-minute timeout. The default seed is 42; pass `--seed 43` to change it.\n\n'
            'Training API (modified local source):\n```python\n'
            'se.set_seed(42, deterministic=True)\n'
            'emb, model = se.train_SpaceExpress_multi(\n'
            '    pair, paths, random_seed=42, deterministic=True, save_model=True)\n```'),
    ]
    nb.metadata.kernelspec = dict(display_name='Python 3', language='python', name='python3')
    path = root / 'seed_reproducibility.ipynb'
    nbformat.write(nb, path)
    client = NotebookClient(nb, timeout=120, resources={'metadata': {'path': str(root)}})
    client.km = client.create_kernel_manager()
    client.km.kernel_spec.argv = [sys.executable, '-m', 'ipykernel_launcher', '-f', '{connection_file}']
    try:
        client.execute()
    finally:
        nbformat.write(nb, path)


def main(args):
    started = time.perf_counter()
    deadline = started + args.max_seconds
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    if (root / 'input').exists() or (root / 'summary.json').exists():
        raise ValueError('Choose a new output directory; existing results are never overwritten')
    (root / 'logs').mkdir(exist_ok=True)
    sys.path.insert(0, str(SOURCE))
    import anndata as ad
    import numpy as np
    import SpaceExpress as se
    runs, comparisons = [], []
    for label, name in DATASETS.items():
        source = args.data_root / name / 'results_v015_retry1'
        inputs = root / 'input' / label
        inputs.mkdir(parents=True)
        provenance = []
        for i in range(2):
            path = source / f'prepared_condition{i}.h5ad'
            backed = ad.read_h5ad(path, backed='r')
            try:
                idx = np.sort(np.random.default_rng(2026 + i).choice(backed.n_obs, 512, replace=False))
                a = backed[idx, :].to_memory()
            finally:
                backed.file.close()
            a.write_h5ad(inputs / f'condition{i}.h5ad')
            provenance.append(dict(source=str(path), shape=list(a.shape), observation_ids=a.obs_names.tolist()))
            se.shortest_path(a.obsm['spatial'], str(inputs / f'shortest{i}.pkl'), k=10)
        (inputs / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
        (root / label).mkdir()
        cases = [
            ('legacy_a', 'installed', args.seed, args.seed, False),
            ('legacy_hash_changed', 'installed', args.seed, 777, False),
            ('fixed_a', 'local', args.seed, args.seed, True),
            ('fixed_b', 'local', args.seed, args.seed, True),
            ('hash_changed', 'local', args.seed, 777, False),
            ('seed_changed', 'local', (args.seed + 1) % 2**32, args.seed, False),
        ]
        for name, implementation, seed, hash_seed, dse in cases:
            command = [sys.executable, '-u', str(Path(__file__).resolve()), '--worker',
                       '--input', str(inputs), '--output', str(root / label / name),
                       '--implementation', implementation, '--seed', str(seed),
                       '--epochs', str(args.epochs), '--device', args.device]
            if dse:
                command.append('--dse')
            env = dict(os.environ, PYTHONHASHSEED=str(hash_seed), CUBLAS_WORKSPACE_CONFIG=':4096:8')
            env.pop('PYTHONPATH', None)
            remaining = deadline - time.perf_counter()
            if remaining <= 0:
                raise TimeoutError('Smoke-test deadline reached')
            print(f'RUN {label}/{name}: seed={seed}, hash_seed={hash_seed}', flush=True)
            with (root / 'logs' / f'{label}_{name}.log').open('w') as log:
                subprocess.run(command, env=env, check=True, stdout=log, stderr=subprocess.STDOUT,
                               timeout=min(180, remaining))
            metadata = json.loads((root / label / f'{name}.json').read_text())
            runs.append(dict(dataset=label, run=name, **metadata))
        comparisons.extend([
            compare(root, label, 'legacy_a', 'legacy_hash_changed', 'observe'),
            compare(root, label, 'fixed_a', 'fixed_b', 'same'),
            compare(root, label, 'fixed_a', 'hash_changed', 'same'),
            compare(root, label, 'fixed_a', 'seed_changed', 'different'),
        ])
    passed = all(c['passed'] and c['exact_dse'] is not False for c in comparisons)
    summary = dict(status='PASS' if passed else 'FAIL', elapsed_seconds=time.perf_counter()-started,
                   runs=runs, comparisons=comparisons)
    (root / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    notebook(root)
    print(json.dumps(summary['comparisons'], indent=2), flush=True)
    print(f"{summary['status']}: {root / 'seed_reproducibility.ipynb'}", flush=True)
    if not passed:
        raise AssertionError('Reproducibility checks failed; see summary.json')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--epochs', type=int, default=30)
    parser.add_argument('--device', choices=['cuda', 'cpu'], default='cuda')
    parser.add_argument('--max-seconds', type=int, default=1200)
    parser.add_argument('--data-root', type=Path, default=HERE.parents[1] / 'Public_Data_Pair_Test_OutlierHVG')
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--input', type=Path)
    parser.add_argument('--implementation', choices=['local', 'installed'], default='local')
    parser.add_argument('--dse', action='store_true')
    args = parser.parse_args()
    if not 0 <= args.seed < 2**32 or not 1 <= args.epochs <= 100:
        parser.error('seed must be in [0, 2**32), epochs must be in [1, 100]')
    worker(args) if args.worker else main(args)
