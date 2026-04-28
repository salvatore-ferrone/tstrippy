# TSTRIPPY Development Plan
Date: 2026-04-28

## Overview

This plan consolidates integrator stabilization, basis-expansion documentation, and composite potential workflows into a unified roadmap. The strategy is to develop features alongside their scientific validation—documentation notebooks serve as both teaching tools and executable convergence tests.

## Completed Work

### Phase 0: Basis-Expansion Foundation ✅
- [x] Implemented axisymmetric Legendre polynomial expansion in Fortran (`mathutils.f90`)
- [x] Created convergence validation in `basis_expansion_verification.ipynb`
- [x] Demonstrated spherical harmonics convergence in `legendre_BFE_orbit_convergence.ipynb`
- [x] Identified spherical harmonics limitations for flattened systems (q < 0.3)
- [x] Refactored tstrippy module structure to allow composite (multi-component) basis potentials
- [x] Started `composite_basis_potential.ipynb` documentation
- [x] Updated `io/` module structure (relocated from `Parsers/`)
- [x] Removed `constants.f90`
- [x] Added Fortran aliasing optimization in `mathutils.f90`

## Active Roadmap

### Phase 1: Bessel Function Expansion (Current Focus)

**Goal:** Validate Bessel function approximations for exponential disk profiles and transition from single-component to composite basis workflows.

#### 1.1: Refactor Exponential Disk Bessel Implementation
- [ ] Separate `bessel_disk_eval_component` logic:
  - Rename to `exponential_disk_bessel_eval_component`
  - Document that Bessel functions solve Poisson equation for arbitrary density profiles (not just this one)
  - Create abstract bessel-function interface in `potentials.f90` for future density profiles
  - Keep exponential disk as canonical first example
- [ ] Update integration tests to reflect new naming
- [ ] Add docstrings explaining when Bessel vs Legendre is preferred:
  - Bessel: disk-like (q < 0.3), flattened systems, arbitrary disk density profiles
  - Legendre: spherical/near-spherical mass distributions

#### 1.2: Convergence Validation Notebook (`bessel_functions_expansion.ipynb`)
- [ ] Single exponential disk profile, systematic convergence study
- [ ] Subsection 1: Density profile accuracy
  - Sweep Bessel function order $n_{\max}$
  - Fix ratio of scale length to scale height (e.g., $h/z_0 = 3, 5, 10$)
  - Show how density profile reproduction improves with order
  - Include comparison to analytic density profile
- [ ] Subsection 2: Potential accuracy
  - Plot convergence of potential vs radius and height
  - Finite-difference gradient check: $\mathbf{a} \stackrel{?}{=} -\nabla\Phi$
- [ ] Subsection 3: Orbit convergence
  - Integrate single test particles in Bessel potential
  - Show conservation metrics (energy drift, angular momentum) vs Bessel order
  - Demonstrate time-reversibility as order increases
  - Iterate over 2–3 different $h/z_0$ ratios to show regime-dependent convergence

#### 1.3: Completion Criteria for Phase 1
- [ ] `bessel_functions_expansion.ipynb` executes without error
- [ ] Convergence plots show monotonic improvement with order
- [ ] Orbit metrics (energy, $L_z$, reversibility) meet target thresholds by $n_{\max} \approx 8\text{–}10$
- [ ] Code comments explain Bessel–Legendre trade-off and guide future density profile implementations
- [ ] Fast test suite passes; no regressions in existing potentials

---

### Phase 2: Composite Basis Potential Validation (Follow-On)

**Goal:** Validate multi-component basis expansions using real Galactic model and globular cluster data.

#### 2.1: Composite Model Setup
- [ ] Create or load Ibata2024 Milky Way model configuration
  - Multiple disk components (thin, thick, or thin + Bessel exponential disk)
  - Halo (Legendre expansion or fixed profile)
  - Bulge (if present in model)
- [ ] Verify component registration and dispatch in `integrator.f90`
- [ ] Test that each component can be toggled independently

#### 2.2: Globular Cluster Orbit Convergence Study (`composite_basis_potential.ipynb`)
- [ ] Load sample of globular clusters from Baumgardt catalog
- [ ] For each cluster:
  - Integrate orbit in full Ibata2024 model
  - Sweep $\ell_{\max}$ (halo) and $n_{\max}$ (disk) independently
  - Record final position, energy, angular momentum
- [ ] Convergence metrics:
  - Orbit endpoint mismatch vs reference (full order)
  - Energy drift vs parameter space
  - Histogram of convergence rates across catalog
- [ ] Identify optimal order trade-off for production runs
- [ ] Document method in notebook with cross-references to component-level validation

#### 2.3: Completion Criteria for Phase 2
- [ ] Notebook executes without error on full catalog subset (e.g., N ≥ 20 clusters)
- [ ] Convergence plateaus well-defined for each component
- [ ] Composite model runtimes are feasible for science production
- [ ] Clear guidance on recommended orders for different science goals

---

## Integrator Stabilization (Background)

The following integrator improvements remain priority for robustness but can proceed in parallel:

- [ ] Linear interpolation for host perturber (replace nearest-time sampling)
- [ ] Consistent time handling at integrator entry points
- [ ] Out-of-range time handling policy
- [ ] Full regression test suite (T0–T4 matrix per prior plan)

**Test Matrix** (from prior session):
- T0: API integrity, unit tests (every push, <10s)
- T1: Numerics correctness, determinism, host interpolation (every push, <20s)
- T2: Time reversibility, timestep convergence (every push or nightly, <45s)
- T3: Bar physics, host perturber regression (nightly, 1–3 min)
- T4: Executable notebooks (nightly, 2–15 min)

---

## Documentation and Testing Philosophy

**Two birds, one stone:** Validation notebooks serve dual purpose:
1. **Executable documentation** – readers understand method and implementation via reproducible examples
2. **Scientific validation** – convergence plots, metric stability, and comparison to theory all in one place
3. **Regression prevention** – notebooks run in CI (or nightly) to catch silent degradation

**Lightweight fast tests** remain in `tests/` for quick feedback.  
**Heavy validation and showcase notebooks** live in `docs/source/` and execute on nightly or pre-release CI.

---

## File Organization

### Current Structure (Post-Refactor)
- `tstrippy/src/potentials.f90` – static potentials (Legendre, Bessel, composites)
- `tstrippy/src/mathutils.f90` – Legendre recursion, Bessel utilities, interpolation, Fortran aliasing
- `tstrippy/src/integrator.f90` – dispatch and integration kernels
- `tstrippy/src/hostperturber.f90`, `galacticbar.f90`, `perturbers.f90` – additional physics
- `tstrippy/io/` – parameter loading and data access (relocated from `Parsers/`)
- `docs/source/basis_expansion_verification.ipynb` – Legendre convergence ✅
- `docs/source/legendre_BFE_orbit_convergence.ipynb` – spherical harmonic orbit validation ✅
- `docs/source/bessel_functions_expansion.ipynb` – Bessel convergence (Phase 1.2)
- `docs/source/composite_basis_potential.ipynb` – multi-component validation (Phase 2.2)

---

## Feedback and Recommendations

The plan is well-structured and pragmatic. Here are notes:

### Strengths
1. **Incremental validation:** Single-component convergence before composite ensures you catch bugs early.
2. **Documentation-as-testing:** Dual-purpose notebooks keep science and code coupled, reducing drift.
3. **Realistic data:** Using Ibata2024 + globular clusters gives concrete performance targets and credibility.
4. **Clear separation of concerns:** Exponential disk → abstract Bessel interface → future profiles is good extensibility.

### Watch-Outs
1. **Naming precision:** Ensure docs clearly state "Bessel functions for Poisson equation" vs "Bessel expansion for exponential disk profile" so future contributors understand the scope and don't conflate method with application.
2. **Convergence plateaus:** In Phase 2 with composite models, you may see order-dependent improvement that stalls due to other approximations (e.g., radial grid resolution, integration timestep). Document those coupling effects explicitly.
3. **Table caching strategy:** When composing multiple Bessel and Legendre expansions, precomputed tables can dominate memory. Consider lazy initialization or disk persistence in `io/` if needed.
4. **CI runtime:** Phase 2 on full cluster catalogs may be slow. Plan for nightly-only execution; keep a small smoke-test subset (e.g., 5 clusters) for push-gate CI.

### Suggestions
1. **Add a "method selector" table** early in `composite_basis_potential.ipynb`:
   - When to use spherical harmonics vs Bessel (q threshold, flattening, profile shape)
   - Computational cost and memory trade-off
2. **Cross-reference at notebook boundaries:** Each convergence notebook should link forward to the next stage so readers understand the scientific flow.
3. **Consider a "best-practices" section** in ARCHITECTURE.md on adding new density profiles to the Bessel framework (for your future self and contributors).

### Overall Assessment
**Direction is sound.** The pipeline from component validation → composite integration → production application is pedagogically and scientifically correct. The dual-purpose documentation approach is modern CI/CD best practice. You're on a good trajectory to produce both robust code and credible science.