#!/usr/bin/env bash
set -euo pipefail

# Kept for old submission scripts. The image already contains the bonsai
# Conda environment, so no runtime download or activation is necessary.
BONSAI_PYTHON=${BONSAI_PYTHON:-/opt/conda/envs/bonsai/bin/python}

if [[ ! -x "$BONSAI_PYTHON" ]]; then
  echo "Bonsai Python not found: $BONSAI_PYTHON" >&2
  echo "Rebuild the MaSIF-mimicry image instead of installing Conda at job runtime." >&2
  exit 1
fi

"$BONSAI_PYTHON" -c 'import sys; assert sys.version_info[:2] == (3, 11)'
echo "Bonsai runtime ready: $BONSAI_PYTHON"
