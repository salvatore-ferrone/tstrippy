# Simulator Redesign (Phase 5)

**Goal**: Rebuild simulator.f90 from the ground up to be extensible, robust, and user-friendly.

## Design Principles

1. **Each module manages its own lifecycle** (gravity, host, perturbers, bar)
2. **Simulator orchestrates** all modules
3. **Forces are commutative**: call order doesn't matter
4. **Pointer-based dispatch**: the simulator builds a registry of active force procedures during setup
5. **Phase space lives in the simulator module**: schemes advance module-held state in place
6. **Trajectory retention is memory-bounded**: only save full trajectories for the first `NparticlesSaved`
7. **Snapshot deferral**: Focus on integration first; snapshot/resumption is Phase 6

## Architecture Overview

### Module Pattern

Each physics module (gravity, host, perturbers, galacticbar) has:
- **Initialization**: Setup-specific to that module (varies per module)
  - Gravity: add components → finalize
  - Host: specify orbit, gravitational model, time evolution
  - Perturbers: specify perturber list, time dependence
  - Bar: specify mass model, orientation/kinematics
- **Finalization**: Allocate/precompute tables if needed
- **Force evaluation**: `compute_force(t, x, y, z, ax, ay, az, phi)`
  - Takes: timestamp t, particle positions (x, y, z)
  - Returns: accelerations (ax, ay, az), optentials (optional: phi)
- **Unbinding check** (host only): track two-body escape energy per particle

### Active Force Registry

The simulator should not branch on `host_active`, `bar_active`, etc. inside the hot force loop.

Instead, during setup/finalization it should build a compact registry of active force evaluators:
- gravity force procedure
- host force procedure
- perturber force procedure
- bar force procedure

Only active modules are inserted into the registry. During integration, the total-force routine iterates over that registry and accumulates contributions. This keeps the hot path simple and makes new modules easier to add.

Suggested pattern:
- define an abstract `force_eval_interface`
- define a derived type like `force_provider_t` that stores a procedure pointer
- keep `nactiveforces` and `active_force_providers(:)` in simulator module state
- rebuild this list when simulator/module configuration changes
- keep reusable scratch arrays in simulator module state so the hot loop does not allocate

For side-effect physics that should run during integration but is not part of the force sum
(for example host-driven unbinding checks), use a second callback registry such as
`active_step_callbacks(:)` rather than inserting special cases into the force loop.

### Integration Scheme Pattern

Each integration scheme (leapfrog, forest_ruth, future schemes):
- **Interface**: `scheme_step()`
    - Takes no phase-space arguments
    - Operates on simulator module variables directly: current positions, velocities, time, timestep
    - Calls the simulator force orchestrator internally
- **Procedure pointer**: default to leapfrog, can switch at runtime
- **Pattern**: Mimic `component_handler_t` from gravity module

### Simulator-Owned State

The simulator module owns and persists:
- Initial conditions: `x0, y0, z0, vx0, vy0, vz0`
- Current phase space during integration: `x, y, z, vx, vy, vz`
- Final phase space after integration: `xf, yf, zf, vxf, vyf, vzf`
- Time bookkeeping: `currenttime, dt, ntimesteps, timestamps`
- Memory-bounded trajectory history for a subset of particles

The integration schemes should treat this state as canonical and advance it in place.

### Simulator Orchestration Flow

```
Python API (via f2py)
    ↓
Simulator module
    ├─ setup phase
    │  ├─ Set G
    │  ├─ Setup gravity (add components, finalize)
    │  ├─ Setup host (if active)
    │  ├─ Setup perturbers (if active)
    │  ├─ Setup bar (if active)
    │  ├─ Build active force registry
    │  └─ Select integration scheme
    ├─ Set initial conditions via `setinitialconditions(...)`
    ├─ Set integration parameters via `setintegrationparameters(t0, dt, nsteps)`
    ├─ Configure optional trajectory retention from memory limit
    └─ Run simulation
        └─ Integrator loop
            │ For each timestep:
            │  ├─ Call gravity force
            │  ├─ (Optionally) Call host force + unbinding check
            │  ├─ (Optionally) Call perturbers force
            │  ├─ (Optionally) Call bar force
            │  ├─ Sum forces → total acceleration
            │  ├─ Advance one step via current scheme pointer
            │  └─ (Optionally) record trajectories for saved particles
            └─ Return final state to Python
```

## Module Responsibilities

### Gravity Module
- Already has component handler pattern
- Simulator uses: `evaluategravityforces()`, `evaluategravitypotential()`
- State: tables for each backend, finalized flag

### Host Module
- Store: orbit trajectory (t, x_host, y_host, z_host, v_host)
- Store: structural parameters, gravitational model
- On force eval: interpolate host position at current t, compute force on particles
- Return: force on each particle from host tidal field
- Also: compute two-body escape energy for each particle, flag unbinding

### Perturbers Module
- Store: list of perturbers (stars, satellites, etc.)
- Store: trajectories or analytic models for each
- On force eval: interpolate/compute each perturber's position at current t
- Return: force on each particle from all perturbers
- Pattern: similar to host but multiple bodies

### Galactic Bar Module
- Store: mass model (parameters), orientation angle (time-dependent)
- On force eval: compute bar orientation at current t, evaluate force
- Return: force on each particle from rotating bar potential

## Integration Scheme Pattern

```fortran
MODULE integration_schemes
    ABSTRACT INTERFACE
        SUBROUTINE scheme_step_interface()
        END SUBROUTINE scheme_step_interface
    END INTERFACE

    PROCEDURE(scheme_step_interface), POINTER :: current_scheme => leapfrog_step

CONTAINS

    SUBROUTINE leapfrog_step()
        ! Advances simulator module state by one dt.
    END SUBROUTINE leapfrog_step

    SUBROUTINE forest_ruth_step()
        ! Advances simulator module state by one dt.
    END SUBROUTINE forest_ruth_step

    SUBROUTINE set_scheme(scheme_name)
        CHARACTER*(*), INTENT(IN) :: scheme_name
        SELECT CASE (scheme_name)
            CASE ("leapfrog")
                current_scheme => leapfrog_step
            CASE ("forest_ruth")
                current_scheme => forest_ruth_step
            CASE DEFAULT
                PRINT*, "WARNING: unknown scheme, keeping current scheme", scheme_name
        END SELECT
    END SUBROUTINE set_scheme

END MODULE integration_schemes
```

## Trajectory Retention

Trajectory storage is distinct from restart snapshots.

- The simulator stores full time histories only for the first `NparticlesSaved` particles.
- `NparticlesSaved` is determined from a user memory limit and the number of timesteps.
- Design rule: `NparticlesSaved < MEM_limit / Nsteps`
- Saved arrays should be allocated once before integration.
- The full system still evolves for all particles; only retained history is truncated.

Suggested stored arrays:
- `xtraj(ntimepoints, nparticles_saved)`
- `ytraj(ntimepoints, nparticles_saved)`
- `ztraj(ntimepoints, nparticles_saved)`
- `vxtraj(ntimepoints, nparticles_saved)`
- `vytraj(ntimepoints, nparticles_saved)`
- `vztraj(ntimepoints, nparticles_saved)`

The first saved index should correspond to the initial conditions, and the last saved index to the final state.

## Force Evaluation Orchestration

Simulator provides a composite force subroutine that:

```fortran
SUBROUTINE evaluate_total_force(t, x, y, z, ax, ay, az, phi)
    REAL*8, INTENT(IN) :: t
    REAL*8, DIMENSION(:), INTENT(IN) :: x, y, z
    REAL*8, DIMENSION(:), INTENT(OUT) :: ax, ay, az
    REAL*8, DIMENSION(:), INTENT(OUT), OPTIONAL :: phi

    ! ax_tmp, ay_tmp, az_tmp, phi_tmp are preallocated module scratch buffers.

    ! Initialize accum
    ax = 0.0d0
    ay = 0.0d0
    az = 0.0d0
    IF (PRESENT(phi)) phi = 0.0d0

    INTEGER :: iforce

    DO iforce = 1, nactiveforces
        CALL active_force_providers(iforce)%compute(t, x, y, z, ax_tmp, ay_tmp, az_tmp, phi_tmp)
        ax = ax + ax_tmp
        ay = ay + ay_tmp
        az = az + az_tmp
        IF (PRESENT(phi)) phi = phi + phi_tmp
    END DO

END SUBROUTINE evaluate_total_force
```

Notes:
- gravity should always be present in the registry once finalized
- host-specific unbinding logic should live in a second active-callback registry that runs after the force sum for the current step
- the hot loop should only iterate over active procedure pointers, not inspect inactive-module flags
- the hot loop should not allocate temporary arrays; scratch storage belongs to simulator module state

## Python-Facing API (f2py wrapping)

**Setup phase**:
```python
tstrippy.simulator.set_gravitational_constant(G)
tstrippy.simulator.add_gravity_component("exponentialdisk", [sigma0, hR, hZ])
tstrippy.simulator.finalize_gravity()

tstrippy.simulator.set_host_orbit(t_array, x_host, y_host, z_host, v_host, ...)
tstrippy.simulator.set_host_gravitational_model("hernquist", [M_host, a_host])

tstrippy.simulator.set_integrator_scheme("leapfrog")  # or "forest_ruth"
```

**Configuration phase**:
```python
tstrippy.simulator.setinitialconditions(x_init, y_init, z_init, vx_init, vy_init, vz_init)
tstrippy.simulator.setintegrationparameters(t0=0.0, dt=0.01, nsteps=1000)
```

**Run phase**:
```python
x, y, z, vx, vy, vz = tstrippy.simulator.run_simulation()
```

## Implementation Order

1. **Design integration scheme abstraction** (zero-argument step procedure pointer)
2. **Refactor simulator-owned state** (initial/current/final phase space)
3. **Extract leapfrog to new step pattern**
4. **Port forest_ruth to new step pattern**
5. **Build active force registry** (procedure-pointer list of active modules)
6. **Add trajectory-retention allocation and writes**
7. **Refactor module activation state** (used during setup, not checked in hot loops)
8. **Wrap in simulator lifecycle** (setup → IC → parameters → run)
9. **Test on each module independently**, then in combinations
10. **Defer restart snapshots to Phase 6**

## Benefits of This Design

✅ **Extensibility**: New integrator schemes and force modules plug in cleanly  
✅ **Clarity**: Each module owns its lifecycle; simulator is orchestrator, not monolith  
✅ **Efficiency**: Inactive modules are skipped (pointer-based dispatch)  
✅ **Robustness**: Forces are commutative; order doesn't matter  
✅ **Maintainability**: Separation of concerns (integration ≠ forces ≠ modules)  
✅ **User-friendly**: Python API is clean and minimal  

## Open Questions (Defer to Later Phases)

- Snapshots: Full state vs. phase-space-only? Deferred to Phase 6.
- Exact memory accounting for trajectory retention: finalize once array shapes are fixed.
- Parallel integration: OpenMP on force eval? Out of scope for Phase 5.
- GPU offload: Out of scope; CPU + OpenMP first.

## Hostcluster Rehaul Contract (New)

This section defines the design contract for replacing legacy host-perturber behavior with a new hostcluster module.

### Goals

1. Keep host state update and force evaluation modular and extensible.
2. Support multiple structural backends (analytic and table-backed).
3. Support time-varying structural parameters without backend-specific branching in simulator hot loops.
4. Integrate cleanly with simulator registry-based force orchestration.

### Layered Hostcluster Design

Hostcluster should be implemented as three internal layers:

1. **Kinematics layer**
    - Owns host orbit arrays and interpolation.
    - Returns host position/velocity at current simulator time.

2. **Structure-evolution layer**
    - Owns structural parameter evolution.
    - Returns `params_current(t)` independent of chosen potential backend.

3. **Backend layer**
    - Consumes `params_current` and relative coordinates.
    - Returns force/potential contribution for all particles.
    - Backends can be analytic (`plummer`, `isochrone`) or table-backed (`king`, later).

### Time-Varying Structure Policy

Hostcluster should provide one parameter provider with three modes:

1. **constant**: fixed parameter vector.
2. **law**: named analytic law with law-specific coefficients (for example, double-exponential mass evolution).
3. **table**: user-provided parameter history interpolated in time.

Recommended precedence when multiple are configured:

- `table > law > constant`

This keeps time evolution backend-agnostic and allows adding new models without changing simulator logic.

### Backend Registration Contract

Use handler registration similar to gravity:

1. backend name
2. expected parameter count
3. force evaluator
4. potential evaluator
5. optional init/finalize hook for table builders/loaders

Initial backend set:

1. `plummer` (analytic)
2. `isochrone` (analytic)
3. `king` (registered scaffold; full table path added in later phase)

### Simulator Interface Contract

Simulator should expose wrapper calls for hostcluster setup:

1. `init_hostcluster_kinematics(...)`
2. `set_hostcluster_backend(model_name)`
3. `set_hostcluster_structure_constant(params)`
4. `set_hostcluster_structure_law(law_name, law_params)`
5. `set_hostcluster_structure_table(times, param_table)`
6. `finalize_hostcluster()`

Force orchestration step behavior:

1. `hostcluster.update_state(currenttime)`
2. `hostcluster.force_on_particles(...)`
3. accumulate into total acceleration/potential

No per-step branch explosion in the simulator hot path.

### Bound/Unbound Tracking Contract

For each particle, define specific relative energy:

$$
E_{\mathrm{rel}} = \frac{1}{2}\|\mathbf{v} - \mathbf{v}_{\mathrm{host}}\|^2 + \Phi_{\mathrm{host}}(r, t)
$$

Bound criterion:

$$
E_{\mathrm{rel}} < -\epsilon
$$

Escape-time policy:

1. require `N_confirm` consecutive unbound checks before confirming escape,
2. store first confirmed escape time,
3. keep hysteresis to avoid noise-driven flip-flops near zero energy.

### Error-Handling Policy

For Python-exposed call paths, avoid `STOP`.

Preferred behavior:

1. warning + no-op for invalid state transitions,
2. warning + return for invalid inputs,
3. explicit lifecycle checks (`clear -> init/set -> finalize -> run`).

## Hostcluster Migration Checklist (Phase-by-Phase)

### Phase 0: Contract Lock

1. Freeze hostcluster public API names and signatures.
2. Freeze parameter-evolution precedence (`table > law > constant`).
3. Freeze backend registration format.

### Phase 1: Scaffolding + Build

1. Add hostcluster module skeleton.
2. Wire simulator wrappers and force accumulation path.
3. Keep only one parity backend (`plummer`) at first.
4. Build + smoke test.

### Phase 2: Analytic Backend Expansion

1. Add `isochrone` backend.
2. Validate `nparams` and backend lookup behavior.
3. Add tests for force/potential sanity and regression.

### Phase 3: Time-Varying Structure

1. Implement `constant` mode.
2. Implement `law` mode (double exponential first).
3. Implement `table` mode interpolation.
4. Add equivalence tests where law/table reduce to constant behavior.

### Phase 4: Bound/Unbound Tracking

1. Implement `E_rel` evaluation.
2. Implement hysteresis and first-escape-time recording.
3. Add tests for edge cases near `E_rel ~= 0`.

### Phase 5: King Backend Scaffold

1. Register `king` backend and lifecycle hooks.
2. Add warning no-op evaluator until table path is ready.
3. Implement table build/load path and validate numerics.

### Phase 6: Hardening

1. Remove remaining `STOP` in Python-facing control paths.
2. Add lifecycle misuse tests.
3. Add mixed-physics simulator integration tests (gravity + hostcluster + optional modules).
