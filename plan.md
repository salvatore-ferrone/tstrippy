# TSTRIPPY Development Plan
Date: 2026-05-05

## Overview

The exponential-disk Bessel effort has now crossed the key production threshold:

- The expensive Bessel/Hankel quadrature is paid offline during component setup
- Runtime force evaluation is table-based and fast
- The interpolation layer has been upgraded to a conservative bicubic-Hermite path so forces are derived from one interpolated potential

Current measured runtime benchmark for the table path is on the order of about `2.5e-07 s / step / particle` with one-time setup of about `0.4-0.5 s` per component at the default table resolution.

The next work is no longer about raw speed rescue. It is now about API cleanup, user-facing validation, and composite-science checks across Legendre + Bessel components.

## API Decisions Locked (2026-05-05)

These are now implementation constraints for the refactor.

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
- Provide explicit unit override call (`setunits(...)`) for non-default workflows.
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

- `cleargravity` — full state reset
- `addgravitycomponent(model_name, ...)` — register one component; invalid `model_name` fails immediately
- `finalizegravity` — validates and freezes all models; builds heavy tables only for components that need them (BFE/table path); analytic-only finalize is a flag flip only
- `finalizegravity` is always required even for analytic-only configurations
- After `finalizegravity`, any mutating call is a hard error until `cleargravity`

#### Component identity

- Components are identified only by integer order of addition (0-indexed internally)
- No user labels; no string-keyed components

#### Independent BFE configuration

- `setsphericalbfedefaults(lmax, nr, r_grid)` — spherical-harmonic BFE settings; analytic and cylindrical components ignore this
- `setcylindricalbfedefaults(nr, nz, nk)` — cylindrical/Bessel-table BFE settings; analytic and spherical components ignore this
- Per-component override deferred to a later version

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

- Strict finalize/config errors with informative messages; code-tagged errors are acceptable.
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

#### 3.1b: Independent backend configuration

- Add clearly separated config paths for:
  - spherical-harmonic BFE settings
  - cylindrical/Bessel-table settings
- Default paths should auto-fill backend defaults unless users opt into expert overrides.

#### 3.1c: Unified component registration

- Introduce a unified component-add surface where users can register analytic or BFE-based components in one composite model.
- Keep backend-specific details internal unless users explicitly request advanced control.

#### 3.2: Make finalize the real handoff point

- `finalizeaxisymmetriccompositebasisexpansion()` should do the expensive precompute needed for runtime evaluation
- Finalize should also make the composite basis integrator-ready
- After finalize, calling `setstaticgalaxy("composite_basis", [G])` should no longer be necessary
- After finalize, mutating configuration calls must hard-error until a full `clear`.

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
- [ ] Single run-oriented API consumes module state (`nsteps`, `nparticles`) without redundant run args
- [ ] Trajectory memory-budget behavior and `nskip` controls are implemented with deterministic defaults
- [ ] `writesnapshot` naming replaces `writestream` in the new API surface
- [ ] Composite API can register analytic and BFE components in one coherent flow
- [ ] Spherical-BFE and cylindrical-BFE hyperparameters are independent and non-conflicting
- [ ] Default BFE hyperparameters are available for simple use and overrideable for expert use

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

#### 4.3: Documentation sequencing (not all today)

- [ ] Notebook-first learning path with increasing complexity:
  - one-orbit integration
  - catalog-scale orbit integration
  - stream generation baseline
  - barred-potential catalog workflows
  - barred-potential stream workflows
- [ ] Keep docs practical and copy-adaptable for user scripts

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

#### 5.3: Preset data packaging (release-aligned)

- [ ] Package core potential preset parameters in read-only form
  - pouliasis2017pii (current)
  - pouliasis2017pi, ibata2024, mcmillan2017 (as implemented)
- [ ] Keep bundled presets read-only to users (copy/modify outside package if needed)
- [ ] Tie bundled-data updates to code releases (no separate data semver for now)

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

Documentation style for near-term release:

- notebook-first
- simple-to-advanced progression
- practical reproducible examples over abstract API catalogs

## File Organization Notes For The Next Refactor

- `tstrippy/src/gravity.f90` (target name; currently `potentials.f90`)
  - owns stateful gravity model definitions and finalized disk runtime evaluators
- `tstrippy/src/simulator.f90` (target name; currently `integrator.f90`)
  - owns simulation state, integration lifecycle, and run execution
- `tstrippy/src/mathutils.f90`
  - should hold only reusable numerical kernels needed by the offline builder or interpolation support
- `docs/source/`
  - should separate method-validation notebooks from production-usage notebooks

## Immediate Next Session Checklist

- [ ] Start refactor in small, testable slices with verification after each change
- [ ] Rename module files in staged steps: `integrator.f90` -> `simulator.f90`, `potentials.f90` -> `gravity.f90`
- [ ] Implement strict lifecycle lock (`clear` -> `init/set` -> `finalize` -> `run`)
- [ ] Introduce single run-oriented API with leapfrog default and optional Forest-Ruth selection
- [ ] Implement trajectory memory budget and `trajectory_nskip` controls
- [ ] Rename `writestream` surface API to `writesnapshot`
- [ ] Implement independent configuration paths for spherical-BFE and cylindrical-BFE settings
- [ ] Add unified composite component registration surface (analytic + BFE)
- [ ] Add Bessel orbit conservation sanity harness and thresholds
- [ ] Add mixed composite validation case: Ibata 2024 halo + exponential-disk table
- [ ] Update user documentation with unit conventions and examples