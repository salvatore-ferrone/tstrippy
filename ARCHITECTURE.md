# Architecture

## Purpose

tstrippy is a hybrid Python/Fortran package for orbit integration and tidal stripping in Milky Way-like potentials.

The architecture is currently in an active refactor from a string-dispatch potential/integrator model toward a stateful gravity plus simulator lifecycle model.

## Current Top-Level API

The top-level package exposes:

- `tstrippy.simulator` (compiled Fortran module)
- `tstrippy.gravity` (compiled Fortran module)
- `tstrippy.mathutils` (compiled Fortran module)
- `tstrippy.io` (pure Python)
- `tstrippy.code` (pure Python helpers)

Import wiring is defined in `tstrippy/__init__.py`.

## Fortran Module Layout

Core Fortran files in `tstrippy/src/`:

- `gravity.f90`: static gravity lifecycle and evaluators
   - lifecycle API: `cleargravity`, `setgravityconstant`, `addgravitycomponent`, `finalizegravity`
   - net evaluators: `evaluategravityforces`, `evaluategravitypotential`
   - basis infrastructure: spherical-harmonic (Legendre) and disk-Bessel/table paths
- `simulator.f90`: integration orchestration and runtime state
   - wraps gravity lifecycle for user-facing simulator workflow
   - assembles total force from gravity + optional perturbers/host/bar/nbody
- `mathutils.f90`: shared numerical primitives
   - Legendre utilities
   - Bessel functions
   - interpolation kernels
- `perturbers.f90`, `hostperturber.f90`, `galacticbar.f90`: additional physics subsystems

## Build Outputs

Meson + f2py builds and installs three extension modules:

- `mathutils`
- `gravity`
- `simulator`

See `meson.build` for exact source sets and wrapping commands.

## Runtime Architecture (Refactor State)

### Gravity lifecycle path (preferred)

1. `cleargravity`
2. optional `setgravityconstant`
3. one or more `addgravitycomponent(model_name, params)`
4. `finalizegravity` (or analytic-only auto-finalize on first evaluate)
5. evaluate via `evaluategravityforces` / `evaluategravitypotential`

This path is the primary API direction and owns `GRAVITY_G` internally.

### Simulator path

`simulator.f90` consumes gravity state through wrapped lifecycle calls and combines with other enabled physics for integration routines.

### Composite BFE path (still parallel)

`gravity.f90` also contains a composite BFE manager (`initcompositegravity`, `addcompositesphericalharmonic...`, `addcompositediskbessel...`, `evaluatecompositegravity`).

This currently coexists with lifecycle `addgravitycomponent` and is an active unification target in the plan.

## Known Architectural Tensions (Tracked in plan.md)

- Lifecycle non-analytic spherical-harmonic state is still shared via module-level `BASIS_*` tables, which can contaminate multi-component configurations.
- Bessel internals are being refactored toward force-only and potential-only paths for consistency with the split evaluator architecture.
- Composite and lifecycle control surfaces overlap and need a canonical mixed analytic+BFE path.

## Build and Validation

Preferred environment:

```bash
conda activate tstrippy
./build.sh
pytest tests/ -q
```

Environment-safe fallback when shell activation is uncertain:

```bash
conda run -n tstrippy ./build.sh
conda run -n tstrippy pytest tests/ -q
```

## Near-Term Direction

The active roadmap (see `plan.md`) is:

1. Complete gravity internal consistency for force/potential split paths.
2. Resolve per-component spherical-harmonic table ownership in lifecycle mode.
3. Decide and enforce one canonical mixed analytic+BFE composition interface.
4. Update documentation and notebooks after API stabilization.

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