#!/bin/bash

# filepath: /home-filer/sferrone/tstrippy/build.sh
# Build and install the package

# rm -rf /obs/sferrone/miniconda3/envs/tstrippy/lib/python3.11/site-packages/tstrippy/

# Clean previous builds
rm -rf builddir

# Configure
meson setup builddir

# Build
meson compile -C builddir

# install 
meson install  -C builddir/
# pip install -e . 