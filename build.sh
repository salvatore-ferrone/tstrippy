#!/bin/bash
echo "PREBUILD SCRIPT"
echo "Check to make sure that the meson compiler is using the conda environment."
echo "a common error is the compiler falls back to an incompatible system compiler."
echo "environment: python: $(which python)"
echo "Environment: gfortran: $(which gfortran)"
echo "Environment: f2py: $(which f2py)"
echo "..."
echo ""



rm -rf builddir
meson setup builddir
meson compile -C builddir
meson install -C builddir/
