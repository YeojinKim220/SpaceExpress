import importlib.util
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


spec = importlib.util.spec_from_file_location('public_pair', Path(__file__).parents[1] / 'examples/public_pair.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def test_prepare_slide_counts_provenance_and_no_input_changes(tmp_path):
    rng = np.random.default_rng(2)
    samples = []
    for i in range(2):
        barcodes = [f'b{x}' for x in range(100)]
        genes = [f'g{x}' for x in range(60)]
        counts = pd.DataFrame(rng.poisson(2, (100, 60)), index=barcodes, columns=genes)
        counts.iloc[0] = 0
        counts.index.name = 'barcode'
        dge = tmp_path / f'dge{i}.csv'
        counts.to_csv(dge)
        coords = pd.DataFrame(rng.uniform(size=(100, 2)), index=barcodes, columns=['x', 'y'])
        coords.index.name = 'barcode'
        coords.iloc[1, 0] = np.nan
        locations = tmp_path / f'locations{i}.csv'
        coords.to_csv(locations)
        labels = tmp_path / f'labels{i}.csv'
        pd.DataFrame({'barcode': barcodes, 'max_cell_type': ['1']*100}).to_csv(labels, index=False)
        samples.append({'name': f's{i}', 'condition': str(i), 'dge': str(dge),
                        'locations': str(locations), 'cell_types': str(labels)})
    config = {'type': 'slide', 'samples': samples, 'target_spots': [80, 80], 'seed': 42, 'n_hvg': 20}
    status = {}
    runner.prepare(config, tmp_path, status)
    assert [s['original_observations'] for s in status['samples']] == [100, 100]
    assert [s['qc_observations'] for s in status['samples']] == [98, 98]
    for i in range(2):
        raw = ad.read_h5ad(tmp_path / f'raw_selected_{i}.h5ad')
        prepared = ad.read_h5ad(tmp_path / f'prepared_condition{i}.h5ad')
        assert raw.shape == (80, 60)
        assert prepared.n_obs == 80
        assert 20 <= prepared.n_vars <= 60
        assert prepared.uns['spaceexpress_preprocessing']['hvg_selection'] == 'per_sample_union'
        assert prepared.uns['spaceexpress_preprocessing']['n_top_genes_per_sample'] == 20
        assert list(raw.obs_names) == list(prepared.obs_names)
        assert prepared.uns['spaceexpress_preprocessing']['mean_sd_outliers_removed']
        assert 'b0' not in raw.obs_names and 'b1' not in raw.obs_names
        np.testing.assert_array_equal(raw[:, prepared.var_names].X.toarray(), prepared.layers['counts'].toarray())

    first = ad.read_h5ad(tmp_path / 'prepared_condition0.h5ad')
    second = ad.read_h5ad(tmp_path / 'prepared_condition1.h5ad')
    assert first.var_names.equals(second.var_names)
    np.testing.assert_array_equal(first.var['highly_variable'] | second.var['highly_variable'], True)
