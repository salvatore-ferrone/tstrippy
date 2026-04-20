# Architecture

## Purpose

tstrippy is a hybrid Python/Fortran package for tidal stripping and orbit integration in Galactic potentials.

The repository is organized around a small Python surface API backed by Fortran kernels compiled into Python extension modules with f2py and Meson.

## High-Level Structure

### User-facing Python package

The top-level package exposes four main entry points:

- `tstrippy.integrator`
- `tstrippy.potentials`
- `tstrippy.ergodic`
- `tstrippy.Parsers`

The import surface is defined in `tstrippy/__init__.py`.

### Fortran computational core

The main numerical work lives in `tstrippy/src/`:

- `potentials.f90`: gravitational potential models and force evaluators
- `integrator.f90`: orbit integration driver and dispatch layer
- `mathutils.f90`: numerical helper utilities
- `perturbers.f90`: additional perturber support
- `hostperturber.f90`: moving host potential support
- `galacticbar.f90`: Galactic bar model support

These Fortran sources are compiled into Python modules with f2py.

### Python support code

- `tstrippy/code/ergodic.py`: particle sampling and related Python-side utilities
- `tstrippy/Parsers/`: parameter loaders and data access helpers
- `tstrippy/data/`: YAML parameter files and packaged reference data

## Build And Packaging

### Build system

The project uses Meson plus f2py, configured in `meson.build`.

Two compiled extension modules are produced:

- `potentials`
- `integrator`

The `integrator` module is built from multiple Fortran source files so that it can dispatch across potentials, perturbers, host perturbers, and bar physics.

### Packaging

Packaging metadata is defined in `pyproject.toml`.

Important current constraints:

- Python requirement is `>=3.9,<3.12`
- Build backend is `mesonpy`
- NumPy is required both for runtime and for f2py-based extension building

## Runtime Data Flow

### Import path

At import time, `tstrippy/__init__.py` promotes selected submodules to the top level:

- `tstrippy.lib.integrator.integrator` becomes `tstrippy.integrator`
- `tstrippy.lib.potentials.potentials` becomes `tstrippy.potentials`
- `tstrippy.code.ergodic` becomes `tstrippy.ergodic`
- `tstrippy.Parsers` remains a Python package namespace

### Typical orbit workflow

A typical user workflow is:

1. Load physical parameters from `tstrippy.Parsers`
2. Choose a static Galactic potential by name
3. Configure the Fortran integrator through `tstrippy.integrator`
4. Set initial conditions
5. Run one of the integration routines
6. Read positions, velocities, timestamps, or stream outputs back into Python

### Potential dispatch

The main dispatch point is `setstaticgalaxy` in `integrator.f90`.

This routine maps a string name such as `"pouliasis2017pii"` or `"miyamotonagai"` to a Fortran procedure pointer. That procedure pointer is then used during time integration.

This means that adding a new potential requires changes in at least two places:

1. implement the potential routine in `potentials.f90`
2. register its string name in `integrator.f90`

## Current Module Responsibilities

### `potentials.f90`

This file contains potential evaluation routines of the form:

- input: parameter vector and arrays of positions
- output: accelerations and potential values

Existing models include:

- Hernquist
- Plummer
- Miyamoto-Nagai
- Allen-Santillan halo
- Long-Murali bar
- Composite Pouliasis et al. Milky Way model

This file is the correct home for any new static gravitational potential model.

### `integrator.f90`

This file owns:

- global integration state
- current particle phase-space arrays
- timestamps
- toggles for enabled physics
- procedure-pointer dispatch for the chosen static potential
- integration routines such as leapfrog and velocity Verlet
- host perturber and perturbing-body wiring

This is the central orchestration layer.

### `mathutils.f90`

This file currently contains only a small interpolation helper.

It is the natural place to grow shared numerical building blocks that are not themselves potentials or integrators, for example:

- interpolation utilities
- Legendre polynomial recurrences
- Gauss-Legendre nodes and weights
- basis-expansion helper functions

### `Parsers`

This package loads parameter files and packaged data into Python-friendly structures.

At present it includes:

- Milky Way potential parameter loaders
- Baumgardt globular cluster catalog helpers
- unit/reference frame helpers

This is the right layer for YAML-backed model configuration.

## Extension Pattern: Adding A New Potential

For a new static potential, the standard path is:

1. Add the evaluator to `tstrippy/src/potentials.f90`
2. Add any numerical helper functions needed in `tstrippy/src/mathutils.f90`
3. Register the model name in `setstaticgalaxy` in `tstrippy/src/integrator.f90`
4. Add a Python-side parser in `tstrippy/Parsers/` if the model has named configuration
5. Add data files in `tstrippy/data/` if the model depends on packaged parameters
6. Add tests and notebook examples before merging

## Planned Basis-Expansion Feature

The current planned feature is an axisymmetric basis-function expansion for a flattened density profile.

Proposed architectural split:

- `mathutils.f90`
  - Legendre polynomial recurrence
  - Legendre derivatives
  - Gauss-Legendre quadrature nodes and weights
  - radial interpolation helpers

- `potentials.f90`
  - basis-expansion potential evaluator
  - storage and initialization of radial coefficient tables
  - runtime evaluation of potential and accelerations from interpolated coefficient tables

- `integrator.f90`
  - dispatch registration for the new potential name

- `Parsers/` and `data/`
  - parameter loading for the basis-expansion model
  - optional table persistence if precomputed coefficient tables are stored on disk

### Expected numerical flow for basis expansion

1. Define density `rho(r, mu)` for an axisymmetric model
2. Project onto Legendre modes to obtain `rho_l(r)`
3. Tabulate `rho_l(r_i)` on a logarithmic radial grid
4. Compute and store `Phi_l(r_i)` and radial derivatives
5. Interpolate those tables at runtime for each particle radius
6. Reconstruct `Phi`, `dPhi/dr`, and angular derivatives from Legendre sums
7. Convert to Cartesian accelerations for the integrator

## Repository Layout

### Root

- `README.md`: public-facing overview and installation
- `pyproject.toml`: package metadata
- `meson.build`: extension build definitions
- `build.sh`: convenience build wrapper
- `tests/`: Python tests and validation notebooks
- `docs/`: Sphinx documentation and notebooks

### Package

- `tstrippy/__init__.py`: top-level API surface
- `tstrippy/src/`: Fortran sources
- `tstrippy/code/`: Python helper modules
- `tstrippy/Parsers/`: parameter/data readers
- `tstrippy/data/`: packaged YAML and FITS data
- `tstrippy/lib/`: installed compiled extension location

## Important Current Constraints

- The build assumes `gfortran` is available
- The package is built around f2py-generated extension modules rather than handwritten Python bindings
- The integrator module owns significant global state, so initialization order matters
- New physics features should be added carefully to avoid breaking the current procedural dispatch model

## Recommended Reading Order For Contributors

If you are new to the codebase, read in this order:

1. `README.md`
2. `tstrippy/__init__.py`
3. `meson.build`
4. `tstrippy/src/potentials.f90`
5. `tstrippy/src/integrator.f90`
6. `tstrippy/Parsers/potential_parameters.py`
7. the relevant notebook or test for the feature you are changing

## Maintenance Notes

This file should be updated whenever one of the following changes:

- a new major Fortran module is introduced
- the dispatch mechanism changes
- the build system changes
- a new public subsystem is added
- the recommended extension pattern changes