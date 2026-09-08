#!/bin/bash
set -euo pipefail
PYTHON=${1:?Python executable}
CONFIG=${2:?Config JSON}
OUTPUT=${3:?Output directory}
STAGE=${4:?Stage}
K=${5:-30}
HERE=$(cd "$(dirname "$0")" && pwd)
export PYTHONHASHSEED=42
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMBA_CACHE_DIR="/tmp/se015-numba-${UID}-${SLURM_JOB_ID:-local}"
export MPLCONFIGDIR="/tmp/se015-mpl-${UID}-${SLURM_JOB_ID:-local}"
export XDG_CACHE_HOME="/tmp/se015-xdg-${UID}-${SLURM_JOB_ID:-local}"
export IPYTHONDIR="/tmp/se015-ipython-${UID}-${SLURM_JOB_ID:-local}"
export MPLBACKEND=Agg
mkdir -p "$OUTPUT"
if [[ "$STAGE" == gpu ]]; then
  nvidia-smi
  /usr/bin/time -v -o "$OUTPUT/probe_time.txt" "$PYTHON" -u "$HERE/public_pair.py" \
    --config "$CONFIG" --output "$OUTPUT" --stage probe
  /usr/bin/time -v -o "$OUTPUT/embed_time.txt" "$PYTHON" -u "$HERE/public_pair.py" \
    --config "$CONFIG" --output "$OUTPUT" --stage embed
else
  /usr/bin/time -v -o "$OUTPUT/${STAGE}_${K}_time.txt" "$PYTHON" -u "$HERE/public_pair.py" \
    --config "$CONFIG" --output "$OUTPUT" --stage "$STAGE" --k "$K"
fi
