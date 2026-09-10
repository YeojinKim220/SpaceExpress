"""Random-state controls for repeatable training on a fixed software/device stack."""

import os
import random

import numpy as np
import torch


def set_seed(seed=42, deterministic=False):
    """Seed Python, NumPy and PyTorch; optionally require deterministic kernels.

    Set PYTHONHASHSEED before starting Python, not inside a running interpreter.
    For deterministic CUDA runs also launch with CUBLAS_WORKSPACE_CONFIG=:4096:8.
    These process-wide settings do not guarantee identical results across GPU
    architectures or PyTorch versions. Unsupported deterministic operations raise.
    """
    if not isinstance(seed, (int, np.integer)) or not 0 <= seed < 2**32:
        raise ValueError('seed must be an integer in [0, 2**32)')
    seed = int(seed)
    if deterministic:
        os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.use_deterministic_algorithms(deterministic)
