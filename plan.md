# TSTRIPPY Development Plan
Date: 2026-05-05, Updated 2026-05-21

## Session Status (2026-06-18)

**Agama backend integrated and validated; next focus is flying spheres / perturbers**

✅ **Completed since the last plan update**:
- Implemented the Agama backend path in `gravity` and `simulator`.
- Added generic Agama model loading via inline INI specs and external INI files.
- Verified the Agama gravity path with targeted smoke and pytest coverage.
- Scaled Agama force/potential outputs by the project gravitational constant.

▶ **Immediate next milestone**:
- Build the new `flyingspheres` perturber module as a multi-object force provider with per-object kinematics and structure tables.
- Keep the first implementation general enough for multiple spherical profiles and future close-encounter diagnostics.

## Session Status (2026-05-21)

**Simulator/Hostcluster contract lock: table-first per-parameter override path**

✅ **Completed this session**:
- Aligned simulator-facing hostcluster API toward `configure_*` naming.
- Locked host model lifecycle requirement: model must be configured before parameter overrides.
- Locked per-parameter override semantics with index-based parameter targeting.
- Locked precedence semantics: per-parameter `table > law > constant`.
- Locked update semantics: latest override wins with warning (no hard failure).
- Locked model-reset semantics: reconfiguring model returns to constant baseline until overrides are re-applied.

⚠️ **Current implementation scope**:
- Table override path is prioritized and testable first.
- Law override path is wired for extensibility but can remain non-evaluating until the next slice.

▶ **Immediate next milestone**:
- Add hostcluster contract tests before physics kernels:
  1. finalize fails without kinematics/model,
  2. invalid `param_index` warnings,
  3. table monotonic-time validation,
  4. override replacement behavior,
  5. model reset clears overrides.

## Session Status (2026-05-14)

**Bessel Phase 4: major far-field normalization issue fixed; residual tails still pending**

✅ **Completed this session**:
- Added diagnostic plotting workflow under `tests/diagnostics/` with:
  - far-field directionality curves,
  - inferred-mass curves,
  - masked Poisson residual map.
- Added/kept reusable regression diagnostics in `tests/test_bessel_physics_validation.py`.
- Fixed far-field closure mass normalization in `tstrippy/src/besselbfe.f90`:
  1. corrected table-based mass quadrature in both R and z,
  2. added analytic closure mass for `exponentialdisk` (`M=2*pi*Sigma0*hR^2`).
- Rebuilt and validated:
  - `conda run -n tstrippy ./build.sh`
  - `conda run -n tstrippy pytest tests/test_bessel_physics_validation.py -q`
  - all 9 tests passing.

⚠️ **Current technical debt**:
- Poisson residual-map tails remain elevated (interior masked p95 still nontrivial), indicating remaining truncation/interpolation/solver-quality work.
- Far-field inferred-mass asymptote is now closure-normalized for exponential disk; this is good for runtime stability, but no longer a pure check of raw transform accuracy at large radius.

▶ **Recommended next milestone (agreed pivot)**:
- Freeze Bessel for now and prioritize simulator completion.
- Add orbit-level energy-conservation diagnostics (`E=T+Phi`) as primary physics sanity gate for integrated trajectories.

## Session Status (2026-05-08)

**Simulator Redesign: Phase 5 INCOMING (from bottom-up)**

Milestone: Planning comprehensive simulator rebuild to ensure extensibility and correctness.

Next Phase (2026-05-08 onwards): Simulator redesign from the ground up with the following requirements:
1. Support multiple integration schemes (leapfrog, forest-ruth, extensible for future schemes)
2. Accept user-provided initial conditions cleanly
3. Interface with multiple force modules (gravity, hostperturber, perturbers, galacticbar) correctly
4. Write resumable snapshots on interrupt
5. Be user-friendly for adding new physics modules

Current status: Pre-design phase—conducting "grill-me" review session to finalize architecture before implementation.

Design corrections locked on 2026-05-08:
- Integration schemes are zero-argument steppers that advance simulator module state in place.
- The simulator owns initial, current, and final phase-space arrays as persistent module state.
- Full trajectory retention is a separate feature from restart snapshots.
- Trajectories are stored only for the first `NparticlesSaved`, chosen from a memory limit and `nsteps`.
- The hot-path force evaluation should iterate over a registry of active procedure pointers, not branch on inactive-module flags.
- Initial conditions and integration parameters remain separate configuration calls: `setinitialconditions(...)` and `setintegrationparameters(...)`.
- Restart snapshots remain deferred until after the core integrator rebuild is working.


## Session Status (2026-05-07)

**Bessel Backend Integration: Phase 3 COMPLETE (with known physics issue)**

Milestone: Bessel backend successfully integrated into gravity dispatch system.

✅ **Completed**:
- Bessel backend (BACKEND_BESSEL=3) registered in handler system
- `exponentialdisk` model wired to bessel backend with 3-param signature
- Two critical bugs fixed:
  1. Memory leak: Added `bessel_clear()` call to `cleargravity()` → no crash on second run
  2. Component indexing: Implemented `bessel_slot_for_component()` helper and fixed `finalizegravity()` bessel loop to use internal slot counter (matching SH pattern)
- Integration test validates: component adds, finalizegravity() completes, force/potential evaluation produces deterministic non-NaN values, second run executes without crash
- All 18 smoke tests pass (no regressions)
- Build succeeds cleanly

⚠️ **Known Issue**:
- Physics values produced by bessel evaluator need validation
- Values are deterministic and non-NaN but magnitude/behavior not yet confirmed correct
- Requires detailed convergence checks and comparison with reference solutions

**Next Session (2026-05-08)**:
- Phase 4: Physics validation for bessel backend
- Use notebooks: basis_expansion_verification, legendre_BFE_orbit_convergence, composite_basis_potential
- May require parameter tuning or table resolution adjustments

### Phase 4 Physics Gates For Bessel Backend

The next Bessel step is not more integration wiring. It is numerical validation with explicit pass/fail gates.

Required validation sequence:

1. Force/potential consistency
  - Check that returned accelerations match finite differences of the returned potential away from the origin and table boundaries.
2. Symmetry checks for even disk density
  - Require `Phi(R,z)=Phi(R,-z)`, `a_R(R,z)=a_R(R,-z)`, `a_z(R,z)=-a_z(R,-z)`, and `a_z(R,0)=0`.
3. Far-field normalization
  - For exponential disk, require asymptotic agreement with the total-mass monopole using `M_tot = 2*pi*Sigma0*hR^2`.
4. Resolution convergence
  - Require convergence under `(NR, NZ, NK)` refinement before trusting the table path.
5. Thick-disk reference comparison
  - Compare the backend against a slow direct reference for the full 3D exponential disk, not just internal consistency.

Observed on 2026-05-08 from the new TEMP validation test:
- The far-field monopole check already fails badly for the current Bessel backend.
- For `Sigma0=1`, `hR=4`, `hZ=0.8` at `R=80`, the backend returns `Phi ~ -1.96e-2` whereas the monopole expectation is `Phi ~ -5.40e-6`.
- This points to a normalization / kernel-form issue, not just a small interpolation error.


- Code-side f2py wrapper issue fixed by shortening spherical projector name:
  - `project_axisymmetric_density_generic` -> `project_axisym_density_generic`
- `gravity.f90` import alias updated to the new projector name.
- `TEMP/build.sh` now builds without `only:` (full public wrapping path).
- Validation passed:
  - `bash build.sh` in `TEMP`
  - `conda run -n tstrippy pytest -q test_gravitymini_smoke.py` (4 passed)

Immediate next microsteps (keep one change per step):

1. Wire `ibata2024halo` into SH evaluation branches in `evaluategravityforcecomponents`, `evaluategravityforces`, and `evaluategravitypotential`.
2. Add one reusable TEMP test that evaluates `ibata2024halo` through gravity force API (not just direct density call).
3. Move SH table construction from evaluate-time to finalize-time in one small slice, then re-test.

## Post-Lunch Roadmap (2026-05-07)

Planned execution order after switching to the development branch:

1. Rebuild from a clean base before further refactors.
  - Recreate/refresh build artifacts from scratch.
  - Confirm import and baseline tests are green before any new feature edits.

2. Continue hardening gravity first (TEMP workflow, then port).
  - Reintroduce remaining analytic models and force evaluators into the new gravity path.
  - Keep additive composition/order-independence guarantees as gating checks.
  - Keep per-component-vs-net-force consistency checks as required regressions.

3. Validate extensibility while re-adding models.
  - Add each analytic model in microsteps.
  - For each model: single-component checks, pairwise composition checks, commutativity checks, and component-sum parity checks.
  - Keep these checks in reusable pytest tests (not one-off scripts).

4. After gravity is stable, introduce a separate Bessel module.
  - Extract Bessel-specific state and evaluators into dedicated module/file(s).
  - Preserve gravity public API and composition semantics.
  - Re-run full gravity regression suite after extraction.

5. Rebuild simulator only after gravity/Bessel stabilization.
  - Reconnect simulator to updated gravity interfaces.
  - Validate integrator outputs against gravity-only reference checks.
  - Defer broader simulator redesign until gravity-side correctness is locked.

Gating rule for this roadmap: no step advances unless build + reusable tests pass at that step.

## Next Session Checklist

1. Run environment-safe validation first:
  - `conda run -n tstrippy ./build.sh`
  - `conda run -n tstrippy pytest tests/ -q`
2. Implement per-component lifecycle output APIs for Phase 2f.2:
  - `evaluategravityforcecomponents(N, x, y, z, force_tensor)`
  - `evaluategravitypotentialcomponents(N, x, y, z, phi_tensor)`
3. Add spherical-harmonic-focused tests:
  - force-only and potential-only parity checks
  - multi-component lifecycle test that exposes/fixes shared `BASIS_*` state contamination
4. Implement Phase 2f.3 per-component non-analytic table state for lifecycle components.
5. Decide Phase 2f.4 canonical mixed analytic+BFE composition path (single user-facing path).
6. After API stabilization, start documentation cleanup and notebook/API rename updates.

## Overview

The exponential-disk Bessel effort has now crossed the key production threshold:

- The expensive Bessel/Hankel quadrature is paid offline during component setup
- Runtime force evaluation is table-based and fast
- The interpolation layer has been upgraded to a conservative bicubic-Hermite path so forces are derived from one interpolated potential

Current measured runtime benchmark for the table path is on the order of about `2.5e-07 s / step / particle` with one-time setup of about `0.4-0.5 s` per component at the default table resolution.

The next work is no longer about raw speed rescue. It is now about API cleanup, user-facing validation, and composite-science checks across Legendre + Bessel components.

## API Decisions Locked (2026-05-05)

These are now implementation constraints for the refactor.

## Execution Contract Update (2026-05-06)

These points were reconfirmed and should be treated as settled for the next implementation slices.

### Scope for the next slices

- Immediate priority is API consistency + extendability for analytic and BFE models.
- `density`, `vcirc`, and `vescape` are explicitly deferred.
- Numerical invariants/conservation checks are post-processing concerns for now (not a gating item for this refactor slice).

### Public UX and lifecycle

- User-facing gravity workflow remains:
  1. clear/init module state
  2. add components by `(model_name, params)`
  3. optional family-level BFE hyperparameter override before finalize
  4. finalize (eager table build)
  5. evaluate
- `gravity` and `simulator` should expose equivalent gravity-facing calls, even though f2py module state is independent.
- Mutating calls after finalize should warn + no-op (never `STOP` in Python-exposed flows).

### Units policy

- Default unit convention remains `(Msun, km/s, kpc)` via default `GRAVITY_G`.
- If user overrides `GRAVITY_G`, unit consistency is user responsibility.

### Component semantics

- Component identity is insertion order only (implicit index).
- Duplicate model names are valid and expected.
- Composite/preset models should be compatible with the same add-component UX and additive composition semantics.

### Refactor sequence (approved)

- First extract spherical-harmonic infrastructure into its own source module/file.
- Then extract Bessel infrastructure into its own source module/file.
- After both extractions, test combined fields that mix spherical-harmonic and Bessel-backed components.

### Open edge-policy item (to lock during implementation)

- Table-bound behavior for BFE lookups:
  - `R > Rmax`: warn and clamp/extend policy to be finalized in code.
  - `R < Rmin`: choose between clamp-at-`Rmin` vs local symmetric handling around zero.

### Lifecycle and naming

- Verb map is standardized:
  - `clear*` full reset
  - `init*` initialize/register module or physics component
  - `set*` assign required numeric/config values
  - `finalize*` validate and hard-lock mutable configuration
  - `run*` execute integration
- Strict lifecycle:
  1. `clear`
  2. `init` and `set`
  3. `finalize`
  4. `run`
- After `finalize`, mutating calls are hard errors until `clear` is called.
- Read-only inspection after `finalize` is allowed.

### Integrator API direction

- Replace multiple integration entrypoints with one run-oriented API.
- Source of truth for `nsteps` and `nparticles` is module state from:
  - initial conditions (`setinitialconditions`-style call)
  - integration parameters (`setintegrationparameters`-style call)
- The run entrypoint should not require redundant `nstep`/`NP` arguments.
- Default integration method is leapfrog; optional method switch supports Forest-Ruth.

### Composite API direction (analytic + BFE unification)

- Composite model API should allow users to add any component through a unified surface, regardless of backend:
  - analytic components
  - spherical-harmonic BFE components
  - cylindrical/Bessel-table BFE components
- User-level component addition should be backend-agnostic at first touch, with backend-specific options exposed only when needed.
- Composite cleanup goals:
  - remove awkward backend-specific init burden from standard workflows
  - maintain clear component identity and ordering
  - preserve strict finalize validation before run

### Independent BFE control surfaces

- Spherical-harmonic BFE controls and cylindrical/Bessel BFE controls must be configured independently.
- Disk-only workflows should not require spherical-BFE parameters.
- Spherical-only workflows should not require cylindrical-table parameters.
- Composite workflows may combine both without parameter namespace collisions.

### Trajectories, snapshots, checkpoints

- Final state is always retained.
- Trajectories are always attempted under a memory budget and may subsample particle count.
- Subsample policy: first `N` particles.
- Initial defaults:
  - trajectories enabled
  - `trajectory_budget_mb = 512`
  - `trajectory_nskip = 1`
  - snapshots disabled unless explicitly enabled
- Snapshot naming is standardized: `writesnapshot` (replace `writestream`).
- Snapshot/checkpoint cadence is step-based (deterministic), not wall-time-based.
- Precision policy:
  - checkpoints: `REAL*8` (restart fidelity)
  - analysis snapshots: `REAL*4` (size/performance)
- Resume scope must include full mutable simulation state:
  - current step/time
  - particle phase-space arrays
  - active integration method
  - loaded module parameters/config flags
  - host/perturber time-index state

### Units and module architecture

- Default external unit convention is fixed to `(kpc, km/s, Msun)`.
- G is stored as explicit module state in the gravity module (`GRAVITY_G`).
- Default value: `G = 4.30091727e-6 kpc (km/s)^2 Msun^-1` (named constant, set at module initialization).
- Override: `setgravityconstant(G)` called before `addgravitycomponent`; hard error to change after `finalizegravity`.
- User-facing `addgravitycomponent` params do **not** include G; G is injected from `GRAVITY_G` at evaluation time.
- `simulator.f90` never handles G directly; it calls gravity evaluators that already embed the correct G.
- No unit conversion layer is needed in Python; all Fortran output is already in the chosen unit system.
- `setunits(...)` is deferred; for non-standard units the user calls `setgravityconstant(G_in_their_units)`.
- Promote gravity subsystem to stateful lifecycle architecture (currently `potentials.f90`, target `gravity.f90`).
- Gravity v1 scope remains force + potential; `density/vcirc/vesc` planned later.

### Gravity module API (locked 2026-05-05)

#### Public evaluator names

Four explicit routines — no overloading:

- `evaluategravityforces(N, x, y, z, ax, ay, az)` — net force, Cartesian
- `evaluategravityforcecomponents(N, x, y, z, ax, ay, az_tensor)` — per-component forces, shape `(NP, 3, NCOMP)`
- `evaluategravitypotential(N, x, y, z, phi)` — net potential, Cartesian
- `evaluategravitypotentialcomponents(N, x, y, z, phi_tensor)` — per-component potential, shape `(NP, NCOMP)`

All evaluation routines take arrays. Scalar (single-particle) input is not supported; caller wraps in length-1 array.

#### Lifecycle

- `cleargravity` — full state reset; resets `GRAVITY_G` to default
- `setgravityconstant(G)` — optional override of G before any `addgravitycomponent`; after `finalizegravity`, emit warning and no-op
- `addgravitycomponent(model_name, ...)` — register one component; G comes from module state, not from user params; invalid `model_name` fails immediately
- `finalizegravity` — validates and freezes all models; builds heavy tables only for components that need them (BFE/table path)
- `finalizegravity` is required for BFE/table configurations; analytic-only configurations may auto-finalize once on first evaluation
- After `finalizegravity`, mutating calls emit warnings and no-op until `cleargravity`

#### Component identity

- Components are identified only by integer order of addition (0-indexed internally)
- No user labels; no string-keyed components

#### Independent BFE configuration

- `setsphericalbfedefaults(lmax, nr, r_grid)` — spherical-harmonic BFE settings; analytic and cylindrical components ignore this
- `setcylindricalbfedefaults(nr, nz, nk)` — cylindrical/Bessel-table BFE settings; analytic and spherical components ignore this
- Per-component override is supported in the direction of the API design; defaults remain the simple path

#### Canonical model names for addgravitycomponent

- `plummer`
- `hernquist`
- `miyamotonagai`
- `longmuralibar`
- `allensantillianhalo`
- `pouliasis2017pii`
- `exponential_oblate_halo`
- `ibata2024halo`
- `exponential_disk_bessel`

#### Method-family naming (documentation and API)

- Use explicit family names: `spherical_harmonic` and `disk_bessel`
- Avoid using `axisymmetric` as the primary user-facing family name

#### Naming refactor contract (must-do)

- Public API must not use `axisymmetric*` naming once the rename slice lands.
- `axisymmetric` can remain only in private/internal comments where mathematically useful.
- Public names must communicate solver family directly:
  - `spherical_harmonic` for Legendre/spherical-harmonic expansion
  - `disk_bessel` for cylindrical/Bessel-table expansion
- Rename policy is direct replacement, not long-lived aliasing.

#### Public rename map (locked for implementation)

- `initaxisymmetricbasisexpansion` -> `initsphericalharmonicbasis`
- `clearaxisymmetricbasisexpansion` -> `clearsphericalharmonicbasis`
- `default_init_basis_expansion` -> `defaultinitsphericalharmonicbasis`
- `axisymmetricbasisexpansion_eval` -> `sphericalharmonicbasis_eval`
- `axisymmetricbasisexpansion_eval_component` -> `sphericalharmonicbasis_eval_component`
- `initaxisymmetriccompositebasisexpansion` -> `initcompositegravity`
- `clearaxisymmetriccompositebasisexpansion` -> `clearcompositegravity`
- `finalizeaxisymmetriccompositebasisexpansion` -> `finalizecompositegravity`
- `axisymmetriccompositebasispotential` -> `evaluatecompositegravity`
- `axisymmetriccompositebasispotential_dispatch` -> `evaluatecompositegravity_dispatch`
- `addcompositeexponentialoblate` -> `addcompositesphericalharmonicexponentialoblate`
- `addcompositeibata2024halo` -> `addcompositesphericalharmonicibata2024halo`
- `addcompositebesselcomponent` -> `addcompositediskbesselcomponent`
- `addcompositebesselexponentialdisk` -> `addcompositediskbesselexponentialdisk`

#### State inspector

- `printgravitystate` — prints all module state to stdout; v1 only, no getter-style array returns yet
- Backend details (table sizes, BFE orders, component kinds) are hidden from normal use but visible via `printgravitystate`

#### Coordinate scope

- All public evaluation is Cartesian only in v1
- Cylindrical or spherical coordinates remain internal implementation details where useful

### Bessel-table defaults (initial locked defaults)

These defaults are intended to be robust for v1 and overrideable by advanced users:

- interpolation backend: conservative bicubic Hermite (required invariant)
- radial table size: `nr = 256`
- vertical table size: `nz = 128`
- radial extent: `R in [1e-3*hR, 200*hR]`
- vertical extent: `|z| in [0, 200*hZ]`
- offline quadrature resolution: `nk = 512`

Design rule:
- keep these defaults for simple API paths
- provide explicit override calls for expert workflows

### Error handling and compatibility

- F2PY safety rule: no hard `STOP` in Python-exposed control paths. Python can hang if Fortran aborts.
- Prefer warning + no-op behavior for invalid state/config/model/parameter flows.
- Informative warning messages remain required; code-tagged warnings are acceptable.
- Backward compatibility is not a requirement for this refactor pass.

### Parallelization policy (locked)

- Near-term policy is inter-simulation parallelization, not intra-simulation threading.
- Primary scaling model:
  - scheduler-level fanout across simulations
  - optional Python multiprocessing inside each scheduler task where useful
- OpenMP/MPI are deferred, not a current requirement.
- Runtime/IO implication:
  - each simulation run must support isolated temporary output directories
  - no shared temp-path assumptions

### Data and reproducibility policy (locked direction)

- Fortran remains responsible for fast runtime writes during integration.
- Python remains responsible for post-processing and packaging user-facing outputs.
- Resume/restart is a first-class goal:
  - checkpoints must be sufficient to resume full mutable simulation state
  - checkpoints are deterministic and step-cadenced
- Output metadata should capture all state needed to reproduce a run.

### File/module naming direction (start point)

- Refactor start point:
  - `integrator.f90` -> `simulator.f90`
  - `potentials.f90` -> `gravity.f90`
- The first implementation slices should keep changes incremental and test after each small step.

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

### 5. Gravity non-analytic component state is currently shared, not per-component

- `evaluategravityforces` and `evaluategravitypotential` loop over components, but Legendre density models (`exponential_oblate_halo`, `ibata2024halo`) rely on module-global `BASIS_*` storage.
- After first initialization, later Legendre components can reuse stale projected tables, so per-component parameters can bleed across component evaluations.
- This blocks correct mixed multi-component use of more than one Legendre density component under `addgravitycomponent`.

### 6. Composite and lifecycle APIs are still parallel control surfaces

- `addgravitycomponent` supports analytic + Legendre-density + direct disk-Bessel paths through `GRAVITY_KIND`.
- `initcompositegravity` / `evaluatecompositegravity` manage an independent `COMPOSITE_*` state with no direct binding into `GRAVITY_KIND` lifecycle dispatch.
- This creates ambiguity for mixed analytic + BFE composites and requires an explicit unification decision.

### 7. Naming-consistency cleanup (resolved)

- [x] `default_init_basis_expansion` renamed to `defaultinitsphericalharmonicbasis`.
- [x] Call sites and private declarations synchronized.

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

### Phase 2b: Gravity Module Lifecycle Refactor ✅
- [x] Promoted `gravity.f90` to the canonical stateful static-gravity module
- [x] Implemented lifecycle entry points: `cleargravity`, `setgravityconstant`, `addgravitycomponent`, `finalizegravity`
- [x] Split public gravity evaluation into explicit force and potential paths:
  - `evaluategravityforces`
  - `evaluategravitypotential`
- [x] Removed legacy combined analytic entrypoints as the public direction for gravity usage
- [x] Established module-owned `GRAVITY_G` so user-facing component params no longer carry `G`
- [x] Verified build success after the gravity refactor

### Phase 2c: Gravity Test Harness ✅
- [x] Added focused lifecycle tests in `tests/test_gravity.py`
- [x] Covered clear/reset, G override, valid/invalid component registration, finalize behavior, force evaluation, and potential evaluation
- [x] Verified the focused gravity suite passes (`24 passed`)
- [x] Updated package-structure expectations away from removed combined gravity entrypoints

### Phase 2d: Unified Gravity Component Surface ✅
- [x] Extended `addgravitycomponent` model support to the locked canonical list:
  - `plummer`, `hernquist`, `miyamotonagai`, `longmuralibar`
  - `allensantillianhalo`, `pouliasis2017pii`
  - `exponential_oblate_halo`, `ibata2024halo`, `exponential_disk_bessel`
- [x] Routed `evaluategravityforces` and `evaluategravitypotential` through the expanded model-kind dispatch
- [x] Implemented analytic-only auto-finalize on first evaluate call
- [x] Kept non-analytic pre-finalize behavior as warning + no-op (no hard `STOP`)
- [x] Updated gravity tests to reflect the agreed lifecycle semantics
- [x] Verified full suite passes (`30 passed`)

### Phase 2e: Naming Refactor And Family Separation (In Progress)

**Goal:** remove ambiguous `axisymmetric` public naming, align public API to `spherical_harmonic` and `disk_bessel`, and keep the state-driven component flow clear.

#### 2e.1: Public naming migration

- [x] Apply the full rename map in `gravity.f90`
- [x] Remove `axisymmetric*` names from public exports and Python-visible wrappers
- [x] Keep only family-explicit public names (`spherical_harmonic*`, `disk_bessel*`, `compositegravity*`)

#### 2e.2: Family-state separation cleanup

- [x] Split naming and comments so spherical-harmonic state and disk-bessel state are clearly distinct
- [ ] Ensure no spherical-harmonic terminology appears in disk-bessel setup/evaluation entrypoints
- [ ] Ensure no disk-bessel terminology appears in spherical-harmonic setup/evaluation entrypoints

#### 2e.3: Component API consistency

- [x] Keep user entrypoint centered on `addgravitycomponent` with canonical model names
- [x] Keep direct evaluator direction consistent with family naming
- [x] Preserve warning + no-op error mode (no hard `STOP` in Python-exposed control paths)

#### 2e.4: Test and documentation updates for rename

- [ ] Update tests and package-structure checks to new names
- [ ] Update docs and notebooks to remove `axisymmetric` public API references
- [x] Verify build + full test suite after rename slice

### Phase 2f: Force/Potential API Divorce For Family Evaluators (Next)

**Goal:** enforce force-only and potential-only internal evaluator paths so user-facing lifecycle calls never compute unnecessary quantities.

#### 2f.1: Spherical-harmonic evaluator split

- [x] Replace combined `sphericalharmonicbasis_eval(N, x, y, z, ax, ay, az, phi_out)` with:
  - `sphericalharmonicbasisforce(N, x, y, z, ax, ay, az)`
  - `sphericalharmonicbasispotential(N, x, y, z, phi_out)`
- [x] Replace combined per-component evaluator with:
  - `sphericalharmonicbasisforce_component(...)`
  - `sphericalharmonicbasispotential_component(...)`
- [x] Keep shared interpolation kernels private so force/potential implementations remain numerically consistent.

#### 2f.2: Gravity lifecycle evaluator split completion

- [x] Ensure `evaluategravityforces` never computes potential temporaries for any component kind.
- [x] Ensure `evaluategravitypotential` never computes force temporaries for any component kind.
- [ ] Introduce per-component lifecycle outputs:
  - `evaluategravityforcecomponents(N, x, y, z, force_tensor)`
  - `evaluategravitypotentialcomponents(N, x, y, z, phi_tensor)`

#### 2f.3: Component-state correctness for non-analytic models

- [ ] Move Legendre table state used by lifecycle components from shared `BASIS_*` storage to per-component storage keyed by component index.
- [ ] Keep composite-only `COMPOSITE_*` state separate unless/until full lifecycle/composite unification is explicitly selected.
- [ ] Validate that multiple Legendre density components with different parameters produce distinct contributions.

#### 2f.4: Composite compatibility decision slice (analytic + BFE)

- [ ] Decide and document one canonical mixed-composite path:
  - Option A: unify through `addgravitycomponent` only
  - Option B: keep `initcompositegravity` path and add explicit analytic-component registration there
- [ ] Eliminate ambiguous dual-path behavior for mixed analytic + BFE workflows.
- [ ] Add tests that cover analytic-only, BFE-only, and mixed analytic+BFE composition semantics.

## Active Roadmap

### Phase 3: Simulator Refactor Around Gravity State ✅

**Goal:** Make `simulator.f90` consume the new stateful `gravity.f90` API directly, without reviving model-dispatch wrappers or duplicating gravity state.

#### 3.1: Lock the simulator/gravity ownership boundary

- `gravity` owns the static galaxy field
- `hostperturber`, `perturbers`, and `galacticbar` keep owning their own state
- `simulator` assembles contributions from all physics modules
- `simulator` should not maintain duplicated truth for whether gravity is configured; it should query gravity/module state directly
- `setstaticgalaxy` remains only as an internal convenience wrapper and is not part of the preferred documented API

#### 3.2: Split simulator assembly into force and potential paths

- Keep the name `HIT`, but redefine it to mean force-only assembly
- Add a separate simulator potential-evaluation routine for the cases that actually need potential output
- Intermediate integrator substeps should call only `HIT`
- Potential should be evaluated only for user-requested timestamps, diagnostics, escape-energy bookkeeping, or explicit output paths

#### 3.3: Make simulator expose the gravity lifecycle directly

- Simulator user-facing gravity calls should mirror gravity naming and semantics as closely as possible:
  - `cleargravitycomponents`
  - `setgravityconstant`
  - `addgravitycomponent`
  - `finalizegravity`
- The user experience through `tstrippy.simulator` should be consistent with direct `tstrippy.gravity` usage
- Extend gravity lifecycle support so named component addition covers the full intended model set, including preset-like models such as `pouliasis2017pii`

#### 3.4: Introduce simulator finalization as the run handoff

- Add `finalizesimulator()` as the main simulator commit step
- `finalizesimulator()` should idempotently call `finalizegravity()` when needed
- The intended simulator lifecycle becomes:
  1. clear
  2. set/init/add physics
  3. finalize simulator
  4. run
- Simulator must refuse to run if the static gravity field is not finalized and no alternative valid force configuration exists

#### 3.5: Normalize the integrator entrypoints around shared assembly

- `leapfrogintime`, `leapfrogtofinalpositions`, and `ruthforestintime` should share the same force-assembly contract
- Repeated setup logic should be reduced only after the new `HIT` and potential-evaluation split is stable
- The future run-oriented API should be built on top of the same shared force/potential assembly routines instead of bypassing them

#### 3.6: Acceptance criteria for Phase 3

- [x] `HIT` is force-only
- [x] A separate simulator potential-evaluation path exists
- [x] Simulator does not duplicate gravity finalized state through `GALAXYISSET`
- [x] `finalizesimulator()` exists and safely forwards to gravity finalization
- [x] `setstaticgalaxy` remains as a hidden one-component wrapper only
- [x] Simulator gravity-facing API is intentionally aligned with gravity naming/behavior
- [x] Full intended gravity models can be added through the lifecycle API

---

### Phase 4: Simulator Test Harness And Integration Stabilization

**Goal:** Build the simulator test suite in slices while the API is being rewritten, so regressions are detected before touching all integration paths.

#### 4.1: First simulator tests

- [ ] Add simulator lifecycle tests for:
  - gravity clear/add/finalize forwarding
  - hidden `setstaticgalaxy` wrapper behavior
  - refusal to run when gravity is not finalized
- [ ] Add a direct consistency test showing simulator static-force assembly matches gravity direct evaluation for a simple analytic model

#### 4.2: Assembly-layer tests before full run tests

- [ ] Test `HIT` as force-only assembly with static gravity only
- [ ] Test combined assembly with host, perturbers, nbody, and bar enabled in isolation-friendly slices
- [ ] Add separate potential-evaluation tests only where simulator physics actually needs potential output

#### 4.3: Integrator tests after assembly is stable

- [ ] Add leapfrog smoke tests on the new assembly contract
- [ ] Add Forest-Ruth smoke tests on the same contract
- [ ] Compare short-time trajectory agreement where both integrators should match
- [ ] Audit repeated initialization/output logic across integration entrypoints and consolidate only after tests exist

#### 4.4: Background cleanup items carried forward

- [ ] Linear interpolation for host perturber (replace nearest-time sampling)
- [ ] Consistent time handling at integrator entry points
- [ ] Out-of-range time handling policy
- [ ] Full regression test suite (`T0` to `T4`)

---

### Phase 5: Composite And Bessel Follow-through

**Goal:** Resume the composite/Bessel cleanup after simulator and gravity lifecycle integration is stable.

#### 5.1: Composite API cleanup

- [ ] Remove mandatory Legendre-style init for disk-only usage
- [ ] Keep spherical-BFE and cylindrical/Bessel configuration independent
- [ ] Preserve a unified component-registration surface across analytic + BFE components
- [ ] Make finalize the true handoff for composite runtime readiness

#### 5.2: Documentation and science validation

- [ ] Update `ARCHITECTURE.md` to reflect the gravity/simulator lifecycle architecture
- [ ] Add Bessel-table orbit sanity and conservation checks
- [ ] Add mixed composite validation cases, including Ibata 2024 + exponential-disk combinations
- [ ] Package preset parameter workflows in a user-facing, reproducible form

## Documentation And Testing Philosophy

Validation notebooks still serve two roles:

1. executable documentation
2. scientific regression tests

But performance now needs an equally explicit place in the workflow. From this point onward, every new disk backend should ship with:

- a scientific comparison against the validated reference path
- a benchmark for `seconds / step / particle`

Documentation style for near-term release:

- notebook-first
- simple-to-advanced progression
- practical reproducible examples over abstract API catalogs

## File Organization Notes For The Current Refactor

- `tstrippy/src/gravity.f90`
  - owns stateful static-gravity model definitions and finalized runtime evaluators
- `tstrippy/src/simulator.f90`
  - owns simulation state, integration lifecycle, and assembly across gravity + other physics modules
- `tstrippy/src/mathutils.f90`
  - should hold only reusable numerical kernels needed by the offline builder or interpolation support
- `docs/source/`
  - should separate method-validation notebooks from production-usage notebooks

## Immediate Next Session Checklist

- [ ] Refactor `HIT` to be force-only while adding a separate simulator potential-evaluation routine
- [ ] Remove duplicated simulator gravity-ready state and query gravity/module state directly
- [ ] Add `finalizesimulator()` and make it safely forward to gravity finalization
- [ ] Extend gravity lifecycle registration to the next required named models, starting with `pouliasis2017pii`
- [ ] Add first simulator tests around lifecycle forwarding and static-force assembly consistency
- [ ] Only after those tests pass, simplify repeated integration setup/output logic across leapfrog and Forest-Ruth