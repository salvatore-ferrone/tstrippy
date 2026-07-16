# Building the python package

## Goal
The goal is to create a python package that leverages fortran code for solving hamilton's equations. Launching the simulations and the i/o is done at the python level. This offers the convience of python with the speed of fortran. 

Getting to this point, has been incredibly difficult and requires much more knowledge about packaging, linking, distribution systems, builds, run-time, etc, than I anticipated. 

Here is a bit of a 'log' of all the things that I learned along the way

## Learning journal

I decided to turn my code into a package because the complexity increased and I needed better organization. Additionally, I was developing on a different machine from where the simulations were being launched. Creating a package and integrating it with Git became an obvious choice. To do so, I started with [Bovy's guide](https://pythonpackaging.info/). This was amazing and introduced me to python packaging basics. His introduction on `why` should I want to package and release my code contains all of the arguments I was already convinced by before setting out on this journey.   However, `numpy` told me to use `meson` to compile my code for building a package. So I switched to that. 

What's the spirit of my code? there is the `simulator` which is the grand orchestrator that has different submodules:

- **simulator.f90**: orchestrates the submodules for the gravitational forces and stores the trajectories of our particles. It solves hamilton's equations.
    - **gravity.f90**: handles the gravitational field of the milky way. It has plug-ins to do so
        - **sphericalharmonics.f90**: solving poisson's equations for slightly flattened systems
        - **besselbfe.f90**: intended for solving poisson's equations for very flattened systems. *This does not work*. I abandoned this and then aimed to use *Agama*
        - **agamabackend.f90**: a wrapper that connects my code to `agama`. 
        - **mathutils.f90**
    - **hostcluster.f90**
        - **mathutils.f90**
    - **perturbers.f90**
        - **mathutils.f90**
    - **galacticbar.f90**
        - **mathutils.f90**

This structure is _okay_. I'm recognizing that it would be nicer if I started with C++. The modules are not classes. Therefore, then can only be instantiated once. If it were OOP, it could be easier to structure. It's not that bad though. Also, I learned fortran derived types, which helps alot. For example, in gravity, I created: 

```fortran
    TYPE, PRIVATE :: component_handler_t
        CHARACTER(LEN=32) :: model_name = ""
        INTEGER :: nparams = 0
        INTEGER :: backend = BACKEND_ANALYTIC
        PROCEDURE(force_eval_iface), POINTER, NOPASS :: force_proc => NULL()
        PROCEDURE(potential_eval_iface), POINTER, NOPASS :: potential_proc => NULL()
        PROCEDURE(density_eval_iface), POINTER, NOPASS :: density_proc => NULL()
    END TYPE component_handler_t
```



Things I've learned:
- source distributions
- wheels
- egg
- shared library
- dynamic links
- darwin 
- operating systems
- the build package
- pip install 
- --no-build-isolation 
- site-packages
- `export`
- editable installs 
- Python Package Authority, The Python Software foundation, The Standard Library
- Conda environments 
- ABI 
- meson versus meson-python 
- numpy.distutils 
- setuptools 
- system prefix 
- *symbols*
- clang flags

What is f2py? 

what does it mean to build a package? 

how are dependencies handled? 

what is meson? 

what is the history of python packaging? 

what is my history with `tstrippy`? 



`f2py` can be used to create `shared-library` objects that can be imported into python. Eventually, I upgraded this to create a python package, so it can be imported in any directory on the system and not just where the *.so file lives. This, has proven to be an enourmous challenge. The first challenge was adapting to the new arm architecture in 2022. This is okay now. However, there is more challenges with the packaging field. The first is the depreciation of `distutils`, and the emergence of `meson` for python packaging. this has many challenges... 

