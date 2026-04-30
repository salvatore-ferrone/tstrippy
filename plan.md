# TSTRIPPY Development Plan
Date: 2026-04-30

## Overview

The exponential-disk Bessel effort has now crossed the key production threshold:

- The expensive Bessel/Hankel quadrature is paid offline during component setup
- Runtime force evaluation is table-based and fast
- The interpolation layer has been upgraded to a conservative bicubic-Hermite path so forces are derived from one interpolated potential

Current measured runtime benchmark for the table path is on the order of about `2.5e-07 s / step / particle` with one-time setup of about `0.4-0.5 s` per component at the default table resolution.

The next work is no longer about raw speed rescue. It is now about API cleanup, user-facing validation, and composite-science checks across Legendre + Bessel components.

## Current Findings

These points were verified by reading `potentials.f90`, `integrator.f90`, and by timing the current workflow.

### 1. Runtime bottleneck is in the force evaluator, not the integrator loop

- `axisymmetriccompositebasispotential` calls `exponential_disk_bessel_eval_component` directly for every Bessel component and every `HIT`.
- `exponential_disk_bessel_eval_component` performs a fresh `nk=256` Gauss-Legendre quadrature for every particle evaluation.
- Inside that quadrature, each node evaluates `exp`, `J0`, and `J1`, so the hot loop is doing expensive transcendental work at every timestep.
- The current path therefore scales like `O(N_particle * N_k * N_step)` with a large constant.

### 2. Composite initialization is carrying Legendre-specific baggage into the disk-only path

- `initaxisymmetriccompositebasisexpansion(G, lmax, r_grid, ncomp)` always allocates Legendre-style radial tables.
- A pure disk/Bessel user still has to provide `lmax` and `r_grid`, even though the current disk evaluator does not consume those tables at runtime.
- That API is wrong for the disk-only case. It exposes implementation details the user should not care about.

### 3. `finalizeaxisymmetriccompositebasisexpansion()` is not acting as a true finalize step

- In `potentials.f90`, finalize currently just checks readiness and flips `COMPOSITE_BASIS_FINALIZED = .TRUE.`.
- It does not precompute a runtime-optimal representation for disk components.
- It does not register the composite basis as the active galaxy in the integrator.
- Users therefore must still call `setstaticgalaxy("composite_basis", [G])`, which duplicates state already implied by the finalized composite object.

### 4. Integrator dispatch is still string-routed instead of state-routed

- `HIT` calls `milkywaypotential(milkwayparams, ...)` only if `GALAXYISSET` was established by `setstaticgalaxy`.
- The composite basis state lives in `potentials.f90`, but the integrator will ignore it unless the extra dispatch registration step is also performed.
- That split responsibility is the source of the current redundancy.

## Completed Work

### Phase 0: Basis-Expansion Foundation ✅
- [x] Implemented axisymmetric Legendre polynomial expansion in Fortran (`mathutils.f90`)
- [x] Created convergence validation in `basis_expansion_verification.ipynb`
- [x] Demonstrated spherical harmonics convergence in `legendre_BFE_orbit_convergence.ipynb`
- [x] Identified spherical harmonics limitations for flattened systems (`q < 0.3`)
- [x] Refactored module structure to allow composite axisymmetric potentials
- [x] Started `composite_basis_potential.ipynb` documentation
- [x] Updated `io/` module structure (relocated from `Parsers/`)
- [x] Removed `constants.f90`
- [x] Added Fortran aliasing optimization in `mathutils.f90`

### Phase 1: Scientific Bessel Prototype ✅
- [x] Refactored the exponential-disk Bessel path toward a generic interface
- [x] Established canonical axisymmetric composite naming in `potentials.f90` and `integrator.f90`
- [x] Preserved backward-compatible wrappers where needed
- [x] Reworked the Bessel quadrature to suppress the earlier ringing problem
- [x] Verified that composite-Bessel orbit integration works scientifically
- [x] Identified that the runtime Bessel design is not production-viable on performance grounds

### Phase 2: Production Table Backend + Conservative Interpolation ✅
- [x] Replaced direct runtime quadrature for exponential-disk composite components with precomputed cylindrical tables
- [x] Moved expensive disk precompute into component setup (offline table build)
- [x] Implemented fast runtime lookup path using table interpolation only
- [x] Upgraded interpolation to conservative bicubic-Hermite evaluation (force from one interpolated potential)
- [x] Added mixed-derivative table support for bicubic patches
- [x] Verified compile + smoke benchmarks for the new path

## Active Roadmap

### Phase 3: Simplify The Composite API Around Finalize (Top Priority)

**Goal:** Make the composite workflow shorter, clearer, and less error-prone for real users.

#### 3.1: Remove mandatory Legendre-style init for disk-only usage

- A user configuring only disk/tabulated components should not need to provide `lmax`
- A user configuring only disk/tabulated components should not need to provide a Legendre radial grid
- Disk components should either:
  - auto-initialize their own `(R,z)` table geometry from physical scale lengths, or
  - accept an optional disk-specific grid override API

`initaxisymmetriccompositebasisexpansion()` should remain only for cases that truly need shared Legendre radial infrastructure.

#### 3.2: Make finalize the real handoff point

- `finalizeaxisymmetriccompositebasisexpansion()` should do the expensive precompute needed for runtime evaluation
- Finalize should also make the composite basis integrator-ready
- After finalize, calling `setstaticgalaxy("composite_basis", [G])` should no longer be necessary

The intended user model is:

1. define components
2. finalize composite basis
3. set initial conditions
4. integrate

Not:

1. define components
2. finalize
3. separately tell the integrator what was already finalized
4. integrate

#### 3.3: Integrator-side architectural change

- The integrator must be able to recognize a finalized composite basis as an active static galaxy without the extra string-dispatch step
- `setstaticgalaxy` should remain for analytic potentials and backward compatibility, but composite basis should no longer depend on it
- The source of truth for composite readiness should be the finalized composite state, not duplicated string registration

#### 3.4: Acceptance criteria for Phase 3

- [ ] Disk-only composite setup does not require `lmax` or a Legendre radial grid
- [ ] Finalize is sufficient to make the composite basis integrator-ready
- [ ] Existing analytic potential workflows through `setstaticgalaxy` still work unchanged
- [ ] Backward-compatible wrappers are preserved where practical, but the new preferred API is shorter and clearer

---

### Phase 4: Bessel Documentation + Orbit Sanity Checks

**Goal:** Package the Bessel-table workflow as reproducible documentation and enforce sanity checks for orbit use.

#### 4.1: Orbit sanity diagnostics (Bessel table path)

- [ ] Add a reproducible conservation script/test for representative disk orbits
- [ ] Add a bounded-error criterion for `dE/E` (max, std, and drift slope)
- [ ] Add circular-orbit and eccentric-orbit sanity cases at multiple radii
- [ ] Add a force/potential consistency check on the interpolation patch outputs

#### 4.2: Documentation updates

- [ ] Update `ARCHITECTURE.md` with the production table + Hermite interpolation design
- [ ] Add/update a user-facing notebook demonstrating recommended setup and units
- [ ] Document known unit pitfalls (`G=1` code units vs physical `G`) with examples

---

### Phase 5: Composite Science Validation (Legendre + Bessel, including Ibata 2024)

**Goal:** Validate mixed-component composite models beyond single-component disk checks.

#### 5.1: Mixed-model validation targets

- [ ] Build and validate a composite model combining Ibata 2024 halo + exponential disk table components
- [ ] Compare orbit morphology and force cuts between mixed-component composite and reference expectations
- [ ] Check robustness across flattening, scale radii, and component mass fraction sweeps
- [ ] Confirm conservation diagnostics remain acceptable in mixed Legendre + Bessel runs

#### 5.2: Notebook/test deliverables

- [ ] Add a dedicated notebook or script for Ibata 2024 + disk composite checks
- [ ] Add a lightweight regression test for mixed composite setup and force-evaluation sanity

---

## Integrator Stabilization (Background)

These remain important, but the disk runtime replacement now has priority because current performance blocks practical use.

- [ ] Linear interpolation for host perturber (replace nearest-time sampling)
- [ ] Consistent time handling at integrator entry points
- [ ] Out-of-range time handling policy
- [ ] Full regression test suite (`T0` to `T4`)

## Documentation And Testing Philosophy

Validation notebooks still serve two roles:

1. executable documentation
2. scientific regression tests

But performance now needs an equally explicit place in the workflow. From this point onward, every new disk backend should ship with:

- a scientific comparison against the validated reference path
- a benchmark for `seconds / step / particle`

## File Organization Notes For The Next Refactor

- `tstrippy/src/potentials.f90`
  - owns the composite component definitions and the finalized disk runtime evaluator
- `tstrippy/src/integrator.f90`
  - owns integration state and should stop requiring redundant composite registration
- `tstrippy/src/mathutils.f90`
  - should hold only reusable numerical kernels needed by the offline builder or interpolation support
- `docs/source/`
  - should separate method-validation notebooks from production-usage notebooks

## Immediate Next Session Checklist

- [ ] Start Phase 3 API cleanup: remove redundant disk-only init requirements where possible
- [ ] Define and document the preferred simplified composite setup sequence
- [ ] Add Bessel orbit conservation sanity harness and thresholds
- [ ] Add mixed composite validation case: Ibata 2024 halo + exponential-disk table
- [ ] Update user documentation with unit conventions and examples