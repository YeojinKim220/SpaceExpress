import importlib
import pickle

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from scipy import sparse
from scipy.spatial.distance import cdist

import SpaceExpress as se

prep_module = importlib.import_module('SpaceExpress.preprocessing')


def make_samples(sparse_input):
    rng = np.random.default_rng(21)
    samples = []
    for i, genes in enumerate((['a', 'b', 'c', 'd', 'only0'], ['only1', 'd', 'c', 'b', 'a'])):
        values = rng.uniform(1, 2, (60 + i * 10, len(genes)))
        values[0, genes.index('a')] = 100
        a = ad.AnnData(sparse.csc_matrix(values) if sparse_input else values,
                       var=pd.DataFrame(index=genes))
        a.layers['original'] = a.X.copy()
        a.obsm['spatial'] = rng.uniform(size=(a.n_obs, 2))
        a.uns['sample_id'] = i
        samples.append(a)
    return samples


def dense(x):
    return x.toarray() if sparse.issparse(x) else np.asarray(x)


@pytest.mark.parametrize('sparse_input', [False, True])
def test_common_outliers_then_sample_hvg_union(monkeypatch, sparse_input):
    samples = make_samples(sparse_input)
    originals = [a.copy() for a in samples]
    common = ['a', 'b', 'c', 'd']
    matrices = [dense(a[:, common].X).copy() for a in samples]
    pooled = np.concatenate(matrices)
    threshold = pooled.mean(axis=0) + 4 * pooled.std(axis=0, ddof=1)
    expected = [np.where(x >= threshold, 0, x) for x in matrices]
    calls = []

    def select(a, n_top_genes, flavor):
        i = a.uns['sample_id']
        calls.append(i)
        assert a.var_names.tolist() == common
        np.testing.assert_allclose(dense(a.X), expected[i])
        assert n_top_genes == 2 and flavor == 'seurat'
        a.var['highly_variable'] = a.var_names.isin(['a', 'b'] if i == 0 else ['b', 'c'])

    monkeypatch.setattr(prep_module.sc.pp, 'highly_variable_genes', select)
    result = se.preprocessing(samples, n_top_genes=2)
    assert calls == [0, 1]
    for i, a in enumerate(result):
        assert a.var_names.tolist() == ['a', 'b', 'c']
        assert a.n_obs == samples[i].n_obs
        np.testing.assert_allclose(dense(a.X), expected[i][:, :3])
        np.testing.assert_array_equal(dense(a.layers['original']), dense(originals[i][:, a.var_names].X))
        np.testing.assert_array_equal(dense(samples[i].X), dense(originals[i].X))
        assert 'spaceexpress_preprocessing' not in samples[i].uns
        meta = a.uns['spaceexpress_preprocessing']
        assert meta['selected_hvg_union'] == 3
        np.testing.assert_array_equal(meta['selected_hvg_per_sample'], [2, 2])
    assert result[0].var['highly_variable'].tolist() == [True, True, False]
    assert result[1].var['highly_variable'].tolist() == [False, True, True]


def test_default_1000_matches_independent_scanpy_union():
    rng = np.random.default_rng(42)
    samples = [ad.AnnData(sparse.csr_matrix(np.log1p(rng.poisson(2, (60, 1500))))) for _ in range(2)]
    pooled = np.concatenate([dense(a.X) for a in samples])
    threshold = pooled.mean(axis=0) + 4 * pooled.std(axis=0, ddof=1)
    expected_sets = []
    for a in samples:
        reference = a.copy()
        x = dense(a.X)
        reference.X = sparse.csr_matrix(np.where(x >= threshold, 0, x).astype(np.float32))
        sc.pp.highly_variable_genes(reference, n_top_genes=1000, flavor='seurat')
        expected_sets.append(set(reference.var_names[reference.var['highly_variable']]))
    result = se.preprocessing(samples)
    expected = sorted(set.union(*expected_sets))
    assert len(expected) > 1000
    for a in result:
        assert a.var_names.tolist() == expected
        assert a.uns['spaceexpress_preprocessing']['n_top_genes_per_sample'] == 1000


def test_small_gene_set_and_inclusive_outlier_boundary(monkeypatch):
    a = ad.AnnData(np.array([[0., 1., 0.], [2., 1., 0.]]),
                   var=pd.DataFrame(index=['varying', 'constant', 'zero']))
    # Across both samples, varying has mean 1 and sample SD sqrt(4/3).
    result = se.preprocessing([a, a.copy()], z_threshold=np.sqrt(3) / 2)
    for item in result:
        assert item.var_names.tolist() == ['constant', 'varying', 'zero']
        np.testing.assert_array_equal(dense(item.X), 0)
        assert item.var['highly_variable'].all()


@pytest.mark.parametrize('n_top_genes', [0, -1, 1.5, True])
def test_invalid_hvg_count(n_top_genes):
    with pytest.raises(ValueError, match='positive integer'):
        se.preprocessing(make_samples(False), n_top_genes=n_top_genes)


@pytest.mark.parametrize('threshold', [0, -1, np.nan, np.inf])
def test_invalid_threshold(threshold):
    with pytest.raises(ValueError, match='finite and positive'):
        se.preprocessing(make_samples(False), z_threshold=threshold)


@pytest.mark.parametrize('case', ['duplicate', 'no_common', 'nan', 'negative', 'empty_sample', 'one_sample'])
def test_invalid_inputs(case):
    pair = make_samples(False)
    if case == 'duplicate':
        pair[0].var_names = ['a'] * pair[0].n_vars
    elif case == 'no_common':
        pair[0].var_names = [f'x{i}' for i in range(pair[0].n_vars)]
    elif case == 'nan':
        pair[0].X[0, 0] = np.nan
    elif case == 'negative':
        pair[0].X[0, 0] = -1
    elif case == 'empty_sample':
        pair[0] = pair[0][:0].copy()
    elif case == 'one_sample':
        pair = pair[:1]
    with pytest.raises(ValueError):
        se.preprocessing(pair)


@pytest.mark.parametrize('multi', [False, True])
def test_training_preserves_preprocessed_union_without_reselecting(tmp_path, monkeypatch, multi):
    samples = se.preprocessing(make_samples(True))
    paths = []
    for i, a in enumerate(samples):
        path = tmp_path / f'paths{i}.pkl'
        distances = cdist(a.obsm['spatial'], a.obsm['spatial'])
        with path.open('wb') as handle:
            pickle.dump({j: [row] for j, row in enumerate(distances)}, handle)
        paths.append(str(path))
    originals = [a.copy() for a in samples]

    def fail_reselection(*args, **kwargs):
        raise AssertionError('Training must preserve the preselected union')

    monkeypatch.setattr(sc.pp, 'highly_variable_genes', fail_reselection)
    kwargs = dict(device='cpu', epochs=1, batch_size=4, num_hvg=1, save_model=True)
    if multi:
        samples[1] = samples[1][:, ::-1].copy()
        emb, model = se.train_SpaceExpress_multi(samples, paths, **kwargs)
        assert [e.shape for e in emb] == [(60, 4), (70, 4)]
    else:
        emb, model = se.train_SpaceExpress(samples[0], paths[0], **kwargs)
        assert emb.shape == (60, 4)
    assert model.encoder[0].in_features == 4
    for a, original in zip(samples, originals):
        np.testing.assert_array_equal(dense(a[:, original.var_names].X), dense(original.X))


def test_training_rejects_mixed_or_incomplete_union():
    samples = se.preprocessing(make_samples(False))
    mixed = [samples[0], make_samples(False)[1]]
    with pytest.raises(ValueError, match='mix prepared'):
        se.train_SpaceExpress_multi(mixed, [])
    samples[1] = samples[1][:, :-1].copy()
    with pytest.raises(ValueError, match='same unique union genes'):
        se.train_SpaceExpress_multi(samples, [])
