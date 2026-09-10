import numpy as np
import pandas as pd
import anndata as ad
import pytest
from scipy import sparse

from SpaceExpress import SpaceExpress_DSE, preprocessing
from SpaceExpress.spaceexpress_dse import empirical_null, spline, robjects


@pytest.mark.parametrize('values', [
    [-1, 1, 1], [-1, -1, -1], [np.nan, np.inf, -np.inf],
    [np.nan, -1, 0, 1, 1], [0, 0, 0], [2],
])
def test_invalid_or_degenerate_statistics_are_not_discoveries(values):
    frame = pd.DataFrame([values], columns=[f'g{i}' for i in range(len(values))])
    result = empirical_null(frame)
    pd.testing.assert_index_equal(result.columns, frame.columns)
    np.testing.assert_array_equal(result.to_numpy(), np.ones((1, len(values))))


def test_fdr_upper_bound_preserves_failed_positions():
    robjects.r('fitdist <- function(...) list(estimate=c(shape=2, rate=1))')
    try:
        values = np.r_[np.linspace(0.1, 0.2, 199), -1, 1, np.nan]
        result = empirical_null(pd.DataFrame([values])).to_numpy()[0]
        assert np.isfinite(result).all()
        assert ((result >= 0) & (result <= 1)).all()
        assert result[199] == result[200] == result[201] == 1
    finally:
        robjects.r('rm(fitdist)')


def test_spline_does_not_repeat_mean_sd_filter():
    rng = np.random.default_rng(12)
    n = 200
    expression = np.zeros(n * 2)
    expression[[4, 13, 70, 205, 216, 240, 333]] = [1, 2, 1, 1, 1, 2, 1]
    frame = pd.DataFrame({'embedding': rng.uniform(size=n*2),
                          'gene': expression, 'group': np.repeat([0, 1], n)})
    current, current_pred, _ = spline(frame, 3)
    assert np.isfinite(current[0]) and current[0] >= 0
    assert np.any(current_pred != 0)


def make_pair(marked):
    rng = np.random.default_rng(7)
    pair = []
    for _ in range(2):
        x = rng.normal(size=(40, 3))
        x[:, 0] = 0
        a = ad.AnnData(sparse.csr_matrix(x), var=pd.DataFrame(index=['constant', 'a', 'b']))
        if marked:
            a.uns['spaceexpress_preprocessing'] = {'mean_sd_outliers_removed': True}
        pair.append(a)
    return pair


@pytest.mark.parametrize('marked', [False, True])
def test_dse_failure_metadata(marked):
    pair = make_pair(marked)
    emb = [np.linspace(0, 1, len(a))[:, None] for a in pair]
    fdr, fitted = SpaceExpress_DSE(emb, pair, k=3, n_jobs=1)
    assert fdr.loc[0, 'constant'] == 1
    for a in fitted:
        assert a.varm['DSE-fit-failed'].loc['constant', 0]
        assert a.varm['DSE-statistic'].loc['constant', 0] == -1


def test_dse_does_not_depend_on_outlier_history():
    pair = make_pair(False)
    pair[0].uns['spaceexpress_preprocessing'] = {'mean_sd_outliers_removed': True}
    emb = [np.linspace(0, 1, len(a))[:, None] for a in pair]
    SpaceExpress_DSE(emb, pair, k=3, n_jobs=1)


def test_preprocessing_marker_survives_h5ad_without_mutating_inputs(tmp_path):
    rng = np.random.default_rng(5)
    pair = [ad.AnnData(sparse.csr_matrix(np.log1p(rng.poisson(2, (80, 30))))),
            ad.AnnData(sparse.csr_matrix(np.log1p(rng.poisson(3, (90, 30)))))]
    before = [a.X.copy() for a in pair]
    cleaned = preprocessing(pair, n_top_genes=10)
    assert cleaned[0].var_names.equals(cleaned[1].var_names)
    assert cleaned[0].n_vars >= 10
    for i, a in enumerate(cleaned):
        assert (before[i] != pair[i].X).nnz == 0
        assert 'spaceexpress_preprocessing' not in pair[i].uns
        path = tmp_path / f'{i}.h5ad'
        a.write_h5ad(path)
        assert ad.read_h5ad(path).uns['spaceexpress_preprocessing']['mean_sd_outliers_removed']


@pytest.mark.parametrize('cell_type', [None, 'cell_type'])
def test_multi_replicate_calls_spline_without_outlier_policy(monkeypatch, cell_type):
    import SpaceExpress.spaceexpress_dse as module
    pair = make_pair(True)
    for a in pair:
        a.obs['cell_type'] = 'type_a'
    calls = []

    def fake_spline(frame, k):
        calls.append(k)
        return np.array([1.0]), np.ones(len(frame)), np.zeros(len(frame))

    name = 'spline_multi_rep_ct' if cell_type else 'spline_multi_rep'
    monkeypatch.setattr(module, name, fake_spline)
    emb = [np.linspace(0, 1, len(a))[:, None] for a in pair]
    SpaceExpress_DSE(emb, pair, k=3, n_jobs=1, multi=True, group_id=[0, 1], cell_type=cell_type)
    assert calls == [3, 3, 3]



@pytest.mark.parametrize('genes', [
    ['b', 'a', 'constant'],  # Same gene set, wrong order.
    ['constant', 'a', 'other'],  # Same shape, different genes.
])
@pytest.mark.parametrize('multi', [False, True])
def test_dse_rejects_misaligned_genes_before_mutating_inputs(genes, multi):
    pair = make_pair(True)
    pair[1].var_names = genes
    emb = [np.linspace(0, 1, len(a))[:, None] for a in pair]
    with pytest.raises(ValueError, match='same genes in the same order'):
        SpaceExpress_DSE(emb, pair, k=3, n_jobs=1, multi=multi,
                         group_id=[0, 1] if multi else None)
    assert all('SpaceExpress' not in a.obsm for a in pair)


@pytest.mark.parametrize('sample_index', [0, 1])
def test_dse_rejects_duplicate_gene_names(sample_index):
    pair = make_pair(True)
    pair[sample_index].var_names = ['constant', 'a', 'a']
    with pytest.raises(ValueError, match='duplicate gene names'):
        SpaceExpress_DSE([np.ones((len(a), 1)) for a in pair], pair, k=3, n_jobs=1)
    assert all('SpaceExpress' not in a.obsm for a in pair)


def test_dse_validates_genes_in_every_replicate():
    samples = make_pair(True) + make_pair(True)
    samples[3].var_names = ['b', 'a', 'constant']
    with pytest.raises(ValueError, match='Sample 3 must have the same genes'):
        SpaceExpress_DSE([np.ones((len(a), 1)) for a in samples], samples,
                         k=3, n_jobs=1, multi=True, group_id=[0, 0, 1, 1])


def test_dse_rejects_empty_gene_set():
    pair = [a[:, :0].copy() for a in make_pair(True)]
    with pytest.raises(ValueError, match='at least one gene'):
        SpaceExpress_DSE([np.ones((len(a), 1)) for a in pair], pair, k=3, n_jobs=1)


@pytest.mark.parametrize('count,multi', [(0, False), (1, False), (3, False), (1, True)])
def test_dse_rejects_invalid_sample_count(count, multi):
    samples = (make_pair(True) + make_pair(True))[:count]
    with pytest.raises(ValueError, match='exactly two samples'):
        SpaceExpress_DSE([np.ones((len(a), 1)) for a in samples], samples,
                         k=3, n_jobs=1, multi=multi)
