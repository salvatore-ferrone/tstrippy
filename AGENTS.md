# AGENTS

## Goal
Help AI coding agents be productive in tstrippy with minimal setup, correct build steps, and safe Fortran refactors.

## Read First
Before making changes, read these project docs:

1. [ARCHITECTURE.md](ARCHITECTURE.md)
2. [plan.md](plan.md)
3. [README.md](README.md)

Prefer linking to these docs in PR notes and chat summaries instead of duplicating long explanations.

## Canonical Environment And Build
1. Activate the conda environment: `conda activate tstrippy`
2. Preferred build command: `./build.sh`
3. Manual equivalent:
   - `meson setup builddir`
   - `meson compile -C builddir`
   - `meson install -C builddir`
4. Quick import check:
   - `python -c "import tstrippy; print(tstrippy.integrator, tstrippy.potentials, tstrippy.mathutils)"`

Notes:
- Keep Python within the supported range in [pyproject.toml](pyproject.toml): >=3.9,<3.12.
- Ensure `gfortran` and `f2py` resolve from the active conda env (not system fallback).

## Test And Validation
1. Run fast tests with `pytest` from repo root.
2. For heavy scientific validation, use notebooks in [docs/source](docs/source) (nightly/manual style, not always push-gated).
3. For basis-expansion checks, prioritize:
   - [docs/source/basis_expansion_verification.ipynb](docs/source/basis_expansion_verification.ipynb)
   - [docs/source/legendre_BFE_orbit_convergence.ipynb](docs/source/legendre_BFE_orbit_convergence.ipynb)
   - [docs/source/composite_basis_potential.ipynb](docs/source/composite_basis_potential.ipynb)

## Source Ownership Boundaries (Fortran)
Use these boundaries to avoid mixing concerns:

1. [tstrippy/src/mathutils.f90](tstrippy/src/mathutils.f90): reusable numerical primitives
   - Legendre recursion
   - Bessel utilities
   - interpolation and aliasing helpers
2. [tstrippy/src/potentials.f90](tstrippy/src/potentials.f90): physical potential model evaluators and model-specific logic
3. [tstrippy/src/integrator.f90](tstrippy/src/integrator.f90): model dispatch and time integration orchestration

If a new potential or renamed evaluator is added in potentials, verify dispatch wiring in integrator.

## Active Refactor Focus: Exponential-Disk Bessel
From [plan.md](plan.md), current focus is separating method from application:

1. Treat Bessel machinery as a generic Poisson-solver method.
2. Keep exponential-disk specifics in model-level evaluators.
3. Rename and refactor toward explicit naming like `exponential_disk_bessel_eval_component`.
4. Update dependent tests/docs/notebooks when names or interfaces change.

When editing both mathutils and potentials in one change:
1. Keep low-level kernels in mathutils.
2. Keep density-profile-specific assembly in potentials.
3. Rebuild and re-run at least fast tests after each coherent refactor step.

## Notebook Workflow For Current Phase
Planned validation notebook: [docs/source/bessel_functions_expansion.ipynb](docs/source/bessel_functions_expansion.ipynb)

Execution expectations:
1. Notebook should run top-to-bottom without hidden state assumptions.
2. Include convergence sweeps for Bessel order and profile shape regimes, as defined in [plan.md](plan.md).
3. Include potential-gradient consistency checks and orbit-convergence diagnostics.
4. Keep cell outputs deterministic where possible for CI/nightly reproducibility.

## Practical Guardrails
1. Do not edit generated build artifact trees unless the task is explicitly about build tooling:
   - `build/`, `builddir/`, `temp.*`, `lib.macosx-*`, `src.macosx-*`
2. Prefer minimal diffs in Fortran files and preserve public entry-point names unless refactor requires rename.
3. After Fortran changes, always rebuild before concluding behavior is correct.

## When In Doubt
If behavior conflicts with assumptions here, use project docs as source of truth and update this file with concise corrections.