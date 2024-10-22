# TSTRIPPY
## TIDAL-STRIPING-PYTHON

![Tidal Stripping](logo.png)

Not quite ready yet, but its getting there. Huge thanks to @JoBovy's [python packaging guide](https://pythonpackaging.info/) and GitHub's Co-Pilot for teaching me how to incorporate Fortran dependencies. 

- Does not work with the new Mac processors.
- only tested on Linux systems
- Started in Jan 2024. 
- requires numpy<=1.22.0 since I'm using numpy.distutils. I don't believe any other package can handle Fortran code. Shoulda probably went with C. 
- need to publish on PyPi

NOTES:
The documentation on `tstrippy.readthedocs.io` is going well. I want to be able to compile the code on readthedoc's computer. However, I'm having a difficult time having it compile the fortran submodules properly. The error messages are long and elusive. 


