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

# export FC=$(which gfortran) 
# this explicity set the fortran compiler, 
# I used this when I had fotran as an argument for the project in the meson.build file "project('tstrippy', ['fortran','c'],version: '0.0.1',license: 'MIT', )"
# having fortran set here prompts meson to search for the fortran compiler based on the shell variable
# so by setting it here, I can use the fortran compiler in the meson.build file
# however, this sucks and I want to make sure that it is handeled by python and conda.
# so now, the commands below in the meson file select the fortran compiler based on the conda environment

# Configure
meson setup builddir 
    # --native-file <(echo "[binaries]"; echo "fortran = 'gfortran'")


# Build
meson compile -C builddir
meson install -C builddir/
