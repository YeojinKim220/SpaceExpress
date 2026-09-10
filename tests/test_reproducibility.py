import pickle
import random

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import torch
from scipy.spatial.distance import cdist

import SpaceExpress as se


def test_set_seed_repeats_all_rngs():
    def draw():
        return random.random(), np.random.rand(), torch.rand(4)
    previous = torch.are_deterministic_algorithms_enabled()
    try:
        se.set_seed(42, deterministic=True)
        a = draw()
        se.set_seed(42, deterministic=True)
        b = draw()
        assert a[:2] == b[:2]
        assert torch.equal(a[2], b[2])
        assert torch.are_deterministic_algorithms_enabled()
    finally:
        torch.use_deterministic_algorithms(previous)


@pytest.mark.parametrize('seed', [-1, 2**32, 0.5])
def test_invalid_seed(seed):
    with pytest.raises(ValueError):
        se.set_seed(seed)


@pytest.mark.parametrize('multi', [False, True])
def test_repeated_training_same_input_and_seed(tmp_path, multi):
    rng = np.random.default_rng(12)
    pair, paths = [], []
    for i in range(2):
        a = ad.AnnData(rng.normal(5, 1, (24, 8)).astype(np.float32),
                       var=pd.DataFrame(index=[f'g{j}' for j in range(8)]))
        a.obsm['spatial'] = rng.uniform(size=(24, 2))
        distances = cdist(a.obsm['spatial'], a.obsm['spatial'])
        path = tmp_path / f'shortest{i}.pkl'
        with path.open('wb') as handle:
            pickle.dump({j: [row] for j, row in enumerate(distances)}, handle)
        pair.append(a)
        paths.append(str(path))
    before = [a.X.copy() for a in pair]
    kwargs = dict(device='cpu', epochs=3, random_seed=42, batch_size=4, deterministic=True)
    previous = torch.are_deterministic_algorithms_enabled()
    try:
        if multi:
            first = se.train_SpaceExpress_multi(pair, paths, **kwargs)
            # Reorder the second sample's genes; name-based alignment must remain stable.
            pair[1] = pair[1][:, ::-1].copy()
            second = se.train_SpaceExpress_multi(pair, paths, **kwargs)
            for a, b in zip(first, second):
                np.testing.assert_array_equal(a, b)
        else:
            first = se.train_SpaceExpress(pair[0], paths[0], **kwargs)
            second = se.train_SpaceExpress(pair[0], paths[0], **kwargs)
            np.testing.assert_array_equal(first, second)
        np.testing.assert_array_equal(pair[0].X, before[0])
    finally:
        torch.use_deterministic_algorithms(previous)
