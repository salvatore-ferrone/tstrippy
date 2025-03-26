# TSTRIPPY
## TIDAL-STRIPING-PYTHON

![Tidal Stripping](logo.png)

Not quite ready yet, but its getting there. Huge thanks to @JoBovy's [python packaging guide](https://pythonpackaging.info/) and GitHub's Co-Pilot for teaching me how to incorporate Fortran dependencies. 

- Does not work with the new Mac processors.
- Tested on Linux systems
- Started in Jan 2024. 
- Compatible with Python 3.12 using meson build system
- Cross-platform build via meson


## Installation
```bash
    conda env create -f environment.yml
    meson setup builddir
    meson compile -C builddir/
    meson install  -C builddir/
```

NOTES:
The documentation on `tstrippy.readthedocs.io` is going well. I want to be able to compile the code on readthedoc's computer. However, I'm having a difficult time having it compile the fortran submodules properly. The error messages are long and elusive. 


