# Architecture

## Purpose

tstrippy is a hybrid Python/Fortran package for tidal stripping and orbit integration in Galactic potentials.

The repository is organized around a small Python surface API backed by Fortran kernels compiled into Python extension modules with f2py and Meson.

## High-Level Structure

### User-facing Python package

The top-level package exposes the following entry points:

**Compiled Fortran modules (if built):**
- `tstrippy.integrator` – orbit integration kernels (compiled from `tstrippy/src/integrator.f90`)
- `tstrippy.potentials` – gravitational potential evaluators (compiled from `tstrippy/src/potentials.f90`)
- `tstrippy.mathutils` – numerical utilities (compiled from `tstrippy/src/mathutils.f90`)

**Pure Python modules:**
- `tstrippy.io` – parameter loaders and data access
- `tstrippy.code` – pure Python support code (lazy-loaded submodules):
  - `tstrippy.code.bfe` – basis-function expansion utilities
  - `tstrippy.code.sampling` – particle sampling and initialization
  - `tstrippy.code.orbits` – orbit analysis utilities

The import surface is defined in `tstrippy/__init__.py`.

**Note on development:** During development before the Fortran modules are compiled, the package will import successfully but issue warnings for missing compiled modules. Full functionality is available once built via `python -m pip install -e . --no-build-isolation`.

### Fortran computational core

The main numerical work lives in `tstrippy/src/`:

- `potentials.f90`: gravitational potential models and force evaluators (static, composite basis, perturbers)
- `integrator.f90`: orbit integration driver and dispatch layer
- `mathutils.f90`: numerical helper utilities (Legendre recursion, Bessel functions, interpolation, aliased arrays)
- `perturbers.f90`: additional perturber support
- `hostperturber.f90`: moving host potential support
- `galacticbar.f90`: Galactic bar model support

These Fortran sources are compiled into Python extension modules with f2py and installed into `tstrippy/lib/`.

### Python support code

- `tstrippy/code/bfe.py` – basis-function expansion Python utilities
- `tstrippy/code/sampling.py` – particle sampling and initialization
- `tstrippy/code/orbits.py` – orbit analysis and utilities
- `tstrippy/io/` – parameter loaders and data access helpers
- `tstrippy/data/` – YAML parameter files and packaged reference data

## Building the Package

**Required before first use or after any changes to Fortran source files.**

### Quick Start

```bash
conda activate tstrippy       # Activate the package environment
./build.sh                    # Compile Fortran modules and install package
```

### What Happens During Build

1. Meson reads `meson.build` and prepares the build configuration
2. f2py wraps the Fortran source files (`potentials.f90`, `integrator.f90`, `mathutils.f90`)
3. gfortran compiles the Fortran code into object files
4. Python extension modules are generated and installed into `tstrippy/lib/`
5. The package is installed in editable mode (`pip install -e .`)

### Post-Build Verification

After building, you should be able to import and use the package:

```python
import tstrippy
print(tstrippy.integrator)    # Should not be None
print(tstrippy.potentials)    # Should not be None
print(tstrippy.mathutils)     # Should not be None
```

If the Fortran modules are missing (e.g., `AttributeError` or `ModuleNotFoundError`), the build did not complete successfully. Check build output for errors.

### Troubleshooting

**Issue:** `ModuleNotFoundError: No module named 'tstrippy.lib.integrator'`  
**Solution:** Run `./build.sh` to compile the Fortran modules.

**Issue:** `gfortran: command not found`  
**Solution:** Install gfortran via Homebrew (`brew install gcc`) or your system package manager.

**Issue:** Build fails with cache or version conflicts  
**Solution:** Clean and rebuild:
```bash
rm -rf build builddir
./build.sh
```

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

This file contains shared numerical building blocks used across the computational core:

- **Legendre polynomial recursion** – axisymmetric basis-function utilities ($P_\ell(\mu)$, $dP_\ell/d\mu$, Gauss-Legendre quadrature)
- **Bessel function utilities** – Bessel-based expansions for disk-like potentials
- **Interpolation helpers** – robust radial interpolation for basis-function coefficient tables
- **Fortran aliasing** – optimized pointer aliasing for efficient array indexing

This is the home for all shared numerical primitives that are not themselves potentials, integrators, or problem-specific routines.

### `io/`

This package loads parameter files and packaged data into Python-friendly structures.

At present it includes:

- Milky Way potential parameter loaders
- Baumgardt globular cluster catalog helpers
- unit/reference frame helpers

This is the right layer for YAML-backed model configuration and data access.

## Extension Pattern: Adding A New Potential

For a new static potential, the standard path is:

1. Add the evaluator to `tstrippy/src/potentials.f90`
2. Add any numerical helper functions needed in `tstrippy/src/mathutils.f90`
3. Register the model name in `setstaticgalaxy` in `tstrippy/src/integrator.f90`
4. Add a Python-side loader in `tstrippy/io/` if the model has named configuration
5. Add data files in `tstrippy/data/` if the model depends on packaged parameters
6. Add tests and notebook examples before merging

## Implemented Basis-Expansion Features

### Axisymmetric Legendre Expansion (✅ Complete)

Spherical harmonic basis for potentials with significant spherical symmetry.

**Flow:**
1. Define density profile $\rho(r, \mu)$ where $\mu = \cos\theta$
2. Project onto Legendre modes to obtain $\rho_\ell(r)$ via Gauss-Legendre quadrature
3. Tabulate $\rho_\ell(r_i)$ on logarithmic radial grid
4. Compute and store $\Phi_\ell(r_i)$ and radial derivatives
5. At runtime, interpolate coefficient tables and reconstruct $\Phi$, $\mathbf{a}$ from Legendre sums

**Validation:** See `docs/source/basis_expansion_verification.ipynb` and `legendre_BFE_orbit_convergence.ipynb`

**Architectural home:**
- `mathutils.f90`: Legendre recursion, derivatives, Gauss-Legendre quadrature
- `potentials.f90`: Legendre potential evaluator and coefficient storage
- `integrator.f90`: dispatch registration
- `io/`: parameter loading and table persistence

### Bessel Expansion for Exponential Disk (✅ Complete, Under Refactoring)

Bessel-function basis for solving Poisson equation in disk-like (flattened, $q \lesssim 0.3$) geometries. Currently implemented for exponential disk density profile; method generalizes to arbitrary axisymmetric disk profiles.

**Flow:**
1. Define disk density $\rho_{\text{disk}}(r, z)$ (e.g., exponential profile)
2. Project onto Bessel modes to obtain $\rho_n(r)$ via Hankel transform
3. Tabulate $\rho_n(r_i)$ on logarithmic radial grid
4. Compute and store $\Phi_n(r_i)$ and radial derivatives
5. At runtime, interpolate coefficient tables and reconstruct $\Phi$, $\mathbf{a}$ from Bessel sums

**Current status:**
- Exponential disk implementation in `potentials.f90` (currently named `bessel_disk_eval_component`)
- **In progress:** refactor to decouple method (Bessel solving Poisson) from application (exponential disk)
  - Rename to `exponential_disk_bessel_eval_component`
  - Create abstract Bessel interface for future density profiles
  - Document method vs application separation clearly

**Validation:** In progress with `docs/source/bessel_functions_expansion.ipynb` (Phase 1.2 in plan)

**Architectural home:**
- `mathutils.f90`: Bessel function utilities, interpolation
- `potentials.f90`: Bessel potential evaluators (one per density profile)
- `integrator.f90`: dispatch registration per profile
- `io/`: parameter loading for each profile variant

### Composite Basis Potentials (✅ Structure in place, validation pending)

Multi-component potentials combining multiple basis expansions (Legendre for halo, Bessel for disks, etc.).

**Architectural innovation:**
- `potentials.f90` now supports repeatable component evaluation
- Integrator can call multiple basis routines per timestep and sum contributions
- Each component maintains independent coefficient tables and parameters
- Allows arbitrary nesting of basis methods

**Validation:** In progress with `docs/source/composite_basis_potential.ipynb` (Phase 2.2 in plan)
- Test case: Ibata2024 Milky Way model with Legendre halo + Bessel disk components
- Application: Globular cluster orbit convergence study using Baumgardt catalog

**Architectural home:**
- `potentials.f90`: multi-component dispatch and composition logic
- `mathutils.f90`: shared interpolation and basis utilities
- `integrator.f90`: per-step multi-component summation
- `io/`: composite model configuration and parameter management

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
- `tstrippy/io/`: parameter/data readers (relocated from `Parsers/`)
- `tstrippy/data/`: packaged YAML and FITS data
- `tstrippy/lib/`: installed compiled extension location

## Important Current Constraints and Design Notes

- The build assumes `gfortran` is available
- The package is built around f2py-generated extension modules rather than handwritten Python bindings
- The integrator module owns significant global state, so initialization order matters
- New physics features should be added carefully to avoid breaking the current procedural dispatch model
- Composite potentials require careful namespacing to avoid coefficient-table collisions; use unique prefixes per component
- Fortran aliasing in `mathutils.f90` provides performance optimization but requires careful memory management to avoid unintended coupling between components

## Best Practices: Adding New Basis-Expansion Density Profiles

If you are extending the Bessel framework to support a new density profile (beyond exponential disk):

1. **Method vs. Application separation:** Bessel functions solve Poisson equations generically. Document the density profile $\rho(r, z)$ separately from the method itself.
   - Name the routine: `<profile_name>_bessel_eval_component` (e.g., `exponential_disk_bessel_eval_component`, `sersic_profile_bessel_eval_component`)
   - Add profile-specific docstring explaining the density form and physical motivation

2. **Coefficient table strategy:**
   - Precompute $\rho_n(r)$ and $\Phi_n(r)$ offline or at initialization
   - Store in separate named arrays to avoid collisions in multi-component configurations
   - Use unique prefixes in `io/` loader configuration files

3. **Validation pipeline:**
   - Create a notebook modeled on `bessel_functions_expansion.ipynb`
   - Show density profile accuracy, potential accuracy (finite-difference test), and orbit conservation across order sweeps
   - Document physical regime (q range, applicability) in notebook and docstrings

4. **Integration into `integrator.f90`:**
   - Register new profile name in `setstaticgalaxy` dispatch
   - Verify coefficient tables are initialized before use
   - Test composition with other components (halo, bulge, etc.)

## Recommended Reading Order For Contributors

If you are new to the codebase, read in this order:

1. `README.md`
2. `tstrippy/__init__.py`
3. `meson.build`
4. `tstrippy/src/potentials.f90`
5. `tstrippy/src/integrator.f90`
6. `tstrippy/io/potential_parameters.py` (or equivalent loader module)
7. the relevant notebook or test for the feature you are changing

If you are working on basis-expansion features, also read:
- `docs/source/basis_expansion_verification.ipynb` (Legendre reference)
- `docs/source/bessel_functions_expansion.ipynb` (Bessel reference, if applicable)

## Maintenance Notes

This file should be updated whenever one of the following changes:

- a new major Fortran module is introduced
- the dispatch mechanism changes
- the build system changes
- a new public subsystem is added
- the recommended extension pattern changes