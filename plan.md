# TSTRIPPY Development Plan
Date: 2026-04-29

## Overview

The scientific side of the exponential-disk prototype is now in acceptable shape: the composite axisymmetric workflow runs, the direct Bessel quadrature no longer shows the earlier ringing pathologies, and single-orbit integrations are producing sensible trajectories.

That is not enough for production. The current runtime path is far too slow for real orbit work:

- Measured composite Bessel benchmark: about `4.74e-04 s / step / particle`
- Target production regime: about `1e-08 s / step / particle`
- Required improvement: about `1e4` to `1e5` in hot-path performance

Because of that, the direct runtime Bessel evaluation is now treated as a scientific prototype, not the production implementation. The next phase is to replace it with a tabulated axisymmetric disk solver whose cost is paid at finalize time, not at every force evaluation.

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

## Active Roadmap

### Phase 2: Replace Runtime Bessel With A Tabulated Cylindrical Solver (Top Priority)

**Goal:** Keep the science of the flattened disk model, but move all expensive work out of the integration hot loop.

#### 2.1: Chosen numerical technique

We will replace the current runtime Bessel-force evaluator with:

- A precomputed axisymmetric cylindrical force table on a 2D `(R, z)` grid
- Runtime evaluation by interpolation only (`Phi(R,z)`, `a_R(R,z)`, `a_z(R,z)`)
- One-time Poisson solve during finalize, not during `HIT`

The preferred construction path is:

- Build a disk-component table at finalize time using a cylindrical Poisson solve that is allowed to be expensive offline
- Use a spectral/Hankel-style precompute only during table construction if it remains the cleanest path
- Store the finalized result as tabulated fields, not as a runtime quadrature evaluator

This is the key architectural decision: even if Hankel/Bessel machinery survives inside the offline builder, it must disappear from the runtime force path.

#### 2.2: Runtime representation requirements

- Disk components must evaluate from a table, not from direct quadrature
- Runtime work per particle should be reduced to interpolation and a small amount of geometry
- The evaluator must return `ax`, `ay`, `az`, and `phi` directly from `(R, z)` with no special functions in the hot loop
- The table format must support future arbitrary axisymmetric disk profiles, not only the exponential disk

#### 2.3: Performance target

- Target at least `1e4` faster than the current direct Bessel path on the single-particle benchmark
- Stretch target: within one to two orders of magnitude of the analytic potentials already in tstrippy
- Benchmark gate for acceptance: report `seconds / step / particle` before and after the replacement on the same orbit setup

#### 2.4: Tomorrow's implementation plan

- [ ] Introduce a new production backend for disk-like axisymmetric components based on tabulated cylindrical fields
- [ ] Keep the current direct Bessel evaluator only as a validation/reference path until the replacement is trusted
- [ ] Move all expensive disk precompute work into finalize
- [ ] Ensure the runtime composite evaluator uses interpolation-only disk components
- [ ] Add a focused benchmark script for `seconds / step / particle`

#### 2.5: Acceptance criteria for Phase 2

- [ ] No direct quadrature or per-call Bessel function evaluation remains in the production hot path
- [ ] Single-particle integration speed improves by at least `1e4` relative to the current prototype benchmark
- [ ] Orbit morphology remains scientifically consistent with the validated prototype
- [ ] Potential-gradient consistency checks pass on the finalized table

---

### Phase 3: Simplify The Composite API Around Finalize (Same Priority As Phase 2)

**Goal:** Remove redundant user-facing setup steps and make the composite workflow reflect the true architecture.

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

### Phase 4: Composite Validation After The Replacement Lands

**Goal:** Re-run the science validation after the production backend is in place.

#### 4.1: Scientific checks

- [ ] Compare orbit morphology between the old Bessel prototype and the new tabulated cylindrical backend
- [ ] Compare `phi`, `a_R`, and `a_z` on representative cuts through the disk
- [ ] Re-run force smoothness checks to ensure the interpolation layer does not reintroduce ringing or grid artifacts
- [ ] Re-run time reversibility and conservation diagnostics for representative disk orbits

#### 4.2: Documentation updates

- [ ] Update `ARCHITECTURE.md` to describe the new production disk backend
- [ ] Reframe the old Bessel notebook as method validation/reference, not the production runtime path
- [ ] Add a short user-facing example showing the simplified composite workflow

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

- [ ] Design the finalized disk table format: stored fields, coordinate system, interpolation order
- [ ] Decide whether the offline builder uses Hankel/FFTLog or another cylindrical Poisson solve internally
- [ ] Refactor finalize so it becomes the precompute-and-activate step
- [ ] Remove the mandatory Legendre-style init from the disk-only user path
- [ ] Add a benchmark harness before further scientific polishing