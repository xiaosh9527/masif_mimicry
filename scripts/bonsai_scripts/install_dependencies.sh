#!/usr/bin/env bash
set -euo pipefail

# Dependencies are installed while building the single project image. This
# script remains as a compatibility entry point and validates that image.
MASIF_PYTHON=${MASIF_PYTHON:-/usr/local/bin/masif-python}
BONSAI_PYTHON=${BONSAI_PYTHON:-/opt/conda/envs/bonsai/bin/python}
EVOEF2_DIR=${EVOEF2_DIR:-/opt/EvoEF2}
STRIDE_EXEC=${STRIDE_EXEC:-/opt/stride/stride}
DEEPTMHMM_DIR=${DEEPTMHMM_DIR:-/opt/deepTMHMM/out}

test -x "$MASIF_PYTHON"
test -x "$BONSAI_PYTHON"
test -x "$EVOEF2_DIR/EvoEF2"
test -x "$STRIDE_EXEC"
test -d "$DEEPTMHMM_DIR"

"$MASIF_PYTHON" -c 'import sys; assert sys.version_info[:2] == (3, 6)'
"$BONSAI_PYTHON" -c 'import sys; assert sys.version_info[:2] == (3, 11); import Bio, numpy, pandas, rdkit'

echo "The bundled MaSIF and bonsai runtimes are ready."
