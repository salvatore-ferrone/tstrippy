echo "QUICK BUILD"


#!/bin/bash
set -euo pipefail

ENV_NAME="tstrippy"

echo "FAST BUILD SCRIPT"
echo "Compiling gravitymini directly with f2py (no meson setup in this repo)."
echo "cwd: $(pwd)"

echo "Using conda env: ${ENV_NAME}"
conda run -n "${ENV_NAME}" python -c "import sys; print('python:', sys.executable)"
conda run -n "${ENV_NAME}" which gfortran || true
conda run -n "${ENV_NAME}" python -c "import numpy; print('numpy:', numpy.__version__)"

if [[ ! -f gravity.f90 || ! -f mathutils.f90 || ! -f sphericalharmonicsbfe.f90 || ! -f besselbfe.f90 ]]; then
    echo "ERROR: expected gravity.f90, mathutils.f90, sphericalharmonicsbfe.f90, and besselbfe.f90 in $(pwd)"
    exit 1
fi

rm -f gravitymini*.so gravitymini*.dylib

conda run -n "${ENV_NAME}" python -m numpy.f2py -c \
    simulator.f90 \
    -m tstrippytest

echo "Build complete: tstrippytest module generated in $(pwd)."
