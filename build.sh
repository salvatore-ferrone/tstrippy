#!/bin/bash

# filepath: /home-filer/sferrone/tstrippy/build.sh
# Build and install the package


# Clean previous builds
rm -rf builddir

# Configure
meson setup builddir

# Build
meson compile -C builddir

# install 
meson install  -C builddir/