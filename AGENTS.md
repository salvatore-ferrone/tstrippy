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
   - `python -c "import tstrippy; print(tstrippy.simulator, tstrippy.gravity, tstrippy.mathutils)"`

If shell activation is uncertain (for example after a fresh terminal/session restore), run through conda explicitly:
- Build: `conda run -n tstrippy ./build.sh`
- Tests: `conda run -n tstrippy pytest tests/ -q`
- Quick import check: `conda run -n tstrippy python -c "import tstrippy; print(tstrippy.simulator, tstrippy.gravity, tstrippy.mathutils)"`

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
2. [tstrippy/src/gravity.f90](tstrippy/src/gravity.f90): static gravity lifecycle, model evaluators, and basis infrastructure
3. [tstrippy/src/simulator.f90](tstrippy/src/simulator.f90): integrator orchestration and coupling to gravity/perturbers

If a new gravity component or renamed evaluator is added in gravity, verify simulator wrapper/dispatch wiring.

## Integration Status: Exponential-Disk Bessel

**Phase 3 Complete (2026-05-07)**: Bessel backend integrated into gravity dispatch system.

- Handler registration: `exponentialdisk` model routes to BACKEND_BESSEL
- Lifecycle: `addgravitycomponent("exponentialdisk", [Sigma0, hR, hZ])` → `finalizegravity()` → evaluate
- Multi-component indexing: Implemented `bessel_slot_for_component()` helper (mirrors SH pattern)
- Memory safety: `bessel_clear()` called in `cleargravity()` for proper lifecycle cleanup
- Validation: Integration test passes; no crashes; deterministic outputs; smoke tests green

**Phase 4 In Progress (updated 2026-05-14)**: Physics validation and convergence analysis.

Current status:
- Bessel backend is deterministic and stable in integration tests.
- Far-field closure normalization bug was fixed in `besselbfe.f90`:
   - corrected table mass quadrature,
   - added analytic closure mass for `exponentialdisk` (`M=2*pi*Sigma0*hR^2`).
- Diagnostic suite and plotting scripts now include inferred-mass and residual-map checks.

Current known issue:
- Poisson residual-map tails remain elevated in parts of the interior masked domain; further numerical-quality work is still pending.

Practical next steps:
1. Pause deep Bessel algorithm changes unless a blocking regression appears.
2. Prioritize simulator completion and orbit-level energy-conservation gates.
3. Keep Bessel diagnostics as regression guards while simulator work proceeds.

## Notebook Workflow For Phase 4
Validation notebooks: [docs/source/bessel_functions_expansion.ipynb](docs/source/bessel_functions_expansion.ipynb)

Execution plan:
1. Top-to-bottom convergence sweeps for Bessel order and profile shape regimes
2. Potential-gradient consistency checks (compare analytical and finite-difference force)
3. Orbit integration tests (compare with reference solutions)
4. Document recommended parameter ranges for accuracy vs speed trade-off
5. Keep cell outputs deterministic for CI/nightly reproducibility

## Practical Guardrails
1. Do not edit generated build artifact trees unless the task is explicitly about build tooling:
   - `build/`, `builddir/`, `temp.*`, `lib.macosx-*`, `src.macosx-*`
2. Prefer minimal diffs in Fortran files and preserve public entry-point names unless refactor requires rename.
3. After Fortran changes, always rebuild before concluding behavior is correct.
4. F2PY safety policy: do not use hard `STOP` in Python-exposed Fortran control paths.
   - Reason: `STOP` can terminate the extension call and leave Python waiting for a response.
   - Use warning + no-op (or status return) for invalid state/config/model/parameter flows.

## When In Doubt
If behavior conflicts with assumptions here, use project docs as source of truth and update this file with concise corrections.