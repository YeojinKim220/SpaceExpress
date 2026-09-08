import numpy as np
import pandas as pd
import anndata as ad
import pytest
from scipy import sparse

from SpaceExpress import SpaceExpress_DSE, select_hvg_after_outlier
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


def test_second_mean_sd_filter_can_be_skipped():
    rng = np.random.default_rng(12)
    n = 200
    expression = np.zeros(n * 2)
    expression[[4, 13, 70, 205, 216, 240, 333]] = [1, 2, 1, 1, 1, 2, 1]
    frame = pd.DataFrame({'embedding': rng.uniform(size=n*2),
                          'gene': expression, 'group': np.repeat([0, 1], n)})
    legacy, legacy_pred, _ = spline(frame, 3)
    current, current_pred, _ = spline(frame, 3, remove_mean_sd_outliers=False)
    assert legacy[0] == -1
    assert np.all(legacy_pred == 0)
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


@pytest.mark.parametrize('marked,expected', [(False, True), (True, False)])
def test_dse_history_and_failure_metadata(marked, expected):
    pair = make_pair(marked)
    emb = [np.linspace(0, 1, len(a))[:, None] for a in pair]
    fdr, fitted = SpaceExpress_DSE(emb, pair, k=3, n_jobs=1)
    assert fdr.loc[0, 'constant'] == 1
    for a in fitted:
        assert a.uns['spaceexpress_dse']['remove_mean_sd_outliers'] is expected
        assert a.varm['DSE-fit-failed'].loc['constant', 0]
        assert a.varm['DSE-statistic'].loc['constant', 0] == -1


def test_mixed_histories_require_explicit_policy():
    pair = make_pair(False)
    pair[0].uns['spaceexpress_preprocessing'] = {'mean_sd_outliers_removed': True}
    with pytest.raises(ValueError, match='Mixed preprocessing'):
        SpaceExpress_DSE([np.ones((40, 1))]*2, pair, k=3, n_jobs=1)


def test_preprocessing_marker_survives_h5ad_without_mutating_inputs(tmp_path):
    rng = np.random.default_rng(5)
    pair = [ad.AnnData(sparse.csr_matrix(np.log1p(rng.poisson(2, (80, 30))))),
            ad.AnnData(sparse.csr_matrix(np.log1p(rng.poisson(3, (90, 30)))))]
    before = [a.X.copy() for a in pair]
    cleaned, genes, _ = select_hvg_after_outlier(pair, n_top_genes=10)
    assert len(genes) == 10
    for i, a in enumerate(cleaned):
        assert (before[i] != pair[i].X).nnz == 0
        assert 'spaceexpress_preprocessing' not in pair[i].uns
        path = tmp_path / f'{i}.h5ad'
        a.write_h5ad(path)
        assert ad.read_h5ad(path).uns['spaceexpress_preprocessing']['mean_sd_outliers_removed']


@pytest.mark.parametrize('cell_type', [None, 'cell_type'])
def test_multi_replicate_forwards_preprocessing_policy(monkeypatch, cell_type):
    import SpaceExpress.spaceexpress_dse as module
    pair = make_pair(True)
    for a in pair:
        a.obs['cell_type'] = 'type_a'
    calls = []

    def fake_spline(frame, k, remove_mean_sd_outliers):
        calls.append(remove_mean_sd_outliers)
        return np.array([1.0]), np.ones(len(frame)), np.zeros(len(frame))

    name = 'spline_multi_rep_ct' if cell_type else 'spline_multi_rep'
    monkeypatch.setattr(module, name, fake_spline)
    emb = [np.linspace(0, 1, len(a))[:, None] for a in pair]
    SpaceExpress_DSE(emb, pair, k=3, n_jobs=1, multi=True, group_id=[0, 1], cell_type=cell_type)
    assert calls == [False, False, False]


def test_explicit_policy_overrides_mixed_history():
    pair = make_pair(False)
    pair[0].uns['spaceexpress_preprocessing'] = {'mean_sd_outliers_removed': True}
    emb = [np.linspace(0, 1, len(a))[:, None] for a in pair]
    _, fitted = SpaceExpress_DSE(emb, pair, k=3, n_jobs=1, remove_mean_sd_outliers=False)
    assert all(not a.uns['spaceexpress_dse']['remove_mean_sd_outliers'] for a in fitted)
