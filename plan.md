# TSTRIPPY Stabilization Plan Before New Features
Date: 2026-04-15

## Goal
Make the code operational, well documented, and up to date before adding new science features.

## Key Technical Notes From This Session
- The host perturber currently selects the nearest host timestamp.  
  This creates force discontinuities and can hurt time-reversibility.
- Linear interpolation in time is the right next step for host force evaluation at leapfrog mid-steps.
- There is likely an initialization-time consistency issue in the leapfrog-to-final-positions path that should be checked first.
- Fast tests should stay very small and deterministic.
- Heavy physics validation should run separately (nightly or pre-release).
- Documentation can absolutely serve as executable scientific validation, but only a lightweight subset should run in routine CI.

## Priority Roadmap

### Phase 1: Integrator Stabilization
1. Ensure consistent time handling at integration start in all integrator entry points.
2. Replace nearest-time host sampling with linear interpolation in time for:
   - Host position
   - Host velocity
   - Host mass evaluation time
3. Add out-of-range time behavior policy:
   - Clamp to boundaries, or
   - Fail fast with a clear error
4. Verify consistency between full-orbit and final-position integrators under identical initial conditions.

### Phase 2: Test Architecture Refresh
1. Split tests into fast and slow tiers.
2. Keep fast tests mandatory on every push.
3. Run slow physics and docs execution on nightly or pre-release jobs.
4. Add regression coverage for the host interpolation and bar-plus-leapfrog behavior.

### Phase 3: Documentation Refresh
1. Keep tutorials lightweight and executable.
2. Keep heavier scientific notebooks as reproducibility artifacts.
3. Use docs to communicate reverse-integrability evidence, while preserving tiny fast assertions in tests.

## Test Matrix

| Tier | Perspective | What it validates | Minimal setup | Pass criteria | Runtime target | Cadence | Marker |
|---|---|---|---|---|---|---|---|
| T0 | Software | Package/API integrity | Import checks | No import/type failures | under 10s | every push | core |
| T0 | Software | Parser and units integrity | Instantiate parsers, verify keys/units | Expected keys and unit types exist | under 10s | every push | core |
| T1 | Software | Potential interface contracts | Small direct calls | Correct output shape and finite values | under 10s | every push | core |
| T1 | Numerics | Integrator state correctness | 1 particle, short run | Array sizes correct, finite outputs, monotonic timestamps | under 15s | every push | core |
| T1 | Numerics | Determinism | Fixed seed repeat run | Identical or tight-tolerance equality | under 15s | every push | core |
| T1 | Numerics | Host interpolation correctness | Synthetic host track with known midpoint | Midpoint interpolation error below tolerance | under 20s | every push | core |
| T1 | Numerics | Integrator path consistency | Same ICs, compare orbit vs final-position path | Final state mismatch below tolerance | under 20s | every push | core |
| T2 | Physics | Time reversibility (static) | N=1 and N=5 short runs | Round-trip relative error under threshold | under 30s | every push | physics_fast |
| T2 | Physics | Timestep convergence | dt, dt/2, dt/4 runs | Error decreases with refinement | under 45s | push or nightly | physics_fast |
| T2 | Physics | Energy drift in static potential | 1 particle, longer run | Relative drift bounded | under 45s | nightly | physics_slow |
| T3 | Physics | Rotating bar invariant behavior | Constant pattern-speed bar case | Invariant drift bounded | 1 to 3 min | nightly | physics_slow |
| T3 | Physics | Host perturber regression | Small stream setup | Stable finite outputs and regression match | 1 to 3 min | nightly | physics_slow |
| T4 | Docs | Executable docs sanity | Selected lightweight notebooks | Notebook execution success | 2 to 5 min | nightly | docs_exec |
| T4 | Docs | Scientific reproducibility | Reverse-integrability notebook set | Key metrics/figures reproduced within tolerance | 5 to 15 min | pre-release | docs_repro |

## Initial Thresholds
- Static reversibility:
  - Median relative error below 1e-8
  - Max relative error below 1e-6
- Orbit-path vs final-position-path mismatch:
  - Below 1e-10 absolute in simple baseline cases
- Static-potential energy drift:
  - Fast tests: max relative drift below 1e-5
  - Slow tests: max relative drift below 1e-6
- Synthetic interpolation midpoint error:
  - Near machine precision for linear test tracks

## Practical Policy Decisions
- Yes: keep tiny particle counts (for example N=5) in mandatory CI tests.
- Yes: keep heavier, publication-style demonstrations in documentation notebooks.
- Yes: treat docs as tests, but only execute a curated fast subset in routine CI.

## Next Working Session Scope
Focus: Implement linear interpolation for host perturber.

Proposed implementation checklist:
1. Add bracket-index search around current time.
2. Compute interpolation weight alpha.
3. Interpolate host position and velocity.
4. Evaluate host mass at current integration time.
5. Use interpolated host state in force evaluation.
6. Add tests:
   - Unit interpolation correctness
   - Leapfrog midpoint force continuity regression
   - Time-reversibility improvement check

## Definition of Done for This Sprint
- Core and physics_fast tiers pass.
- Host interpolation regression tests pass.
- No discontinuous force artifacts from nearest-time switching.
- Documentation page updated with a short validation result and method note.