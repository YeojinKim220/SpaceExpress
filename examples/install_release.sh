#!/bin/bash
set -euo pipefail
BASE_PYTHON=${1:?Pass a Python 3.11+ executable with the scientific/R environment}
ENV_DIR=${2:?Pass a new environment directory}
if [[ -e "$ENV_DIR" ]]; then
  echo "Refusing to overwrite an existing environment: $ENV_DIR" >&2
  exit 1
fi
"$BASE_PYTHON" -m venv --system-site-packages "$ENV_DIR"
"$ENV_DIR/bin/python" -m pip install --index-url https://pypi.org/simple \
  --no-cache-dir --report "$ENV_DIR/pip_install_report.json" 'spaceexpress[notebooks]==0.1.5'
"$ENV_DIR/bin/python" -c 'import importlib.metadata as m, SpaceExpress; assert m.version("spaceexpress") == "0.1.5"; print(m.version("spaceexpress"), SpaceExpress.__file__)'
"$ENV_DIR/bin/python" -m pip freeze > "$ENV_DIR/pip_freeze.txt"
"$ENV_DIR/bin/python" -m pip check
