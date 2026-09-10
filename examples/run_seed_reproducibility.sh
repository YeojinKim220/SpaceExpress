#!/bin/bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
SCRATCH=$(cd "$HERE/../../.." && pwd)
BASE_ENV=${BASE_ENV:-/storage/home/hcoda1/9/ykim3030/.conda/envs/spaceexpress-env-test}
PYTHON=${PYTHON:-$SCRATCH/envs/spaceexpress-pypi-0.1.5/bin/python}
OUTPUT="$HERE/../../Reproducibility_Seed_Test"
SEED=42
FOREGROUND=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --output) OUTPUT=$2; shift 2 ;;
    --seed) SEED=$2; shift 2 ;;
    --foreground) FOREGROUND=1; shift ;;
    *) echo "Unknown option: $1" >&2; exit 64 ;;
  esac
done
OUTPUT=$(realpath -m "$OUTPUT")
if [[ -e "$OUTPUT/launch.log" || -e "$OUTPUT/summary.json" || -e "$OUTPUT/input" ]]; then
  echo "Choose a new --output directory: $OUTPUT" >&2
  exit 1
fi
mkdir -p "$OUTPUT"
export PATH="$(dirname "$PYTHON"):$BASE_ENV/bin:$PATH"
export R_HOME="$BASE_ENV/lib/R"
export LD_LIBRARY_PATH="$BASE_ENV/lib:${LD_LIBRARY_PATH:-}"
export PYTHONHASHSEED="$SEED"
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLBACKEND=Agg
export MPLCONFIGDIR="/tmp/se-seed-mpl-$UID-$$"
export NUMBA_CACHE_DIR="/tmp/se-seed-numba-$UID-$$"
export IPYTHONDIR="/tmp/se-seed-ipython-$UID-$$"
export JUPYTER_RUNTIME_DIR="/tmp/se-seed-jupyter-$UID-$$"
if [[ "$FOREGROUND" == 1 ]]; then
  exec timeout --kill-after=10s 25m "$PYTHON" -u "$HERE/seed_reproducibility.py" --output "$OUTPUT" --seed "$SEED"
fi
nohup timeout --kill-after=10s 25m "$PYTHON" -u "$HERE/seed_reproducibility.py" \
  --output "$OUTPUT" --seed "$SEED" > "$OUTPUT/launch.log" 2>&1 < /dev/null &
printf '%s\n' "$!" > "$OUTPUT/launch.pid"
printf 'PID: %s\nLog: %s/launch.log\n' "$!" "$OUTPUT"
