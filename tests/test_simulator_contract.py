import numpy as np
import pytest


tstrippy = pytest.importorskip("tstrippy")
sim = tstrippy.simulator
gravity = tstrippy.gravity

DEFAULT_MB = 1024.0


def _set_basic_problem(n_particles=32, nsteps=100):
    _configure_plummer_gravity()

    x = np.linspace(0.0, 1.0, n_particles)
    y = np.linspace(1.0, 2.0, n_particles)
    z = np.linspace(2.0, 3.0, n_particles)
    vx = np.zeros(n_particles)
    vy = np.zeros(n_particles)
    vz = np.zeros(n_particles)

    sim.setinitialconditions(x, y, z, vx, vy, vz)
    sim.setscheme("leapfrog", [0.0, 1.0e-2, nsteps])


@pytest.fixture(autouse=True)
def _reset_simulator_state():
    sim.clear()
    sim.cleargravitycomponents()
    gravity.clear()
    sim.max_ram_mb = DEFAULT_MB
    sim.orbit_ram_limit_mb = DEFAULT_MB
    yield
    sim.clear()
    sim.cleargravitycomponents()
    gravity.clear()
    sim.max_ram_mb = DEFAULT_MB
    sim.orbit_ram_limit_mb = DEFAULT_MB


def _configure_plummer_gravity(mass=1.0, scale_radius=1.0, g=1.0):
    sim.cleargravitycomponents()
    sim.set_gravitational_constant(g)
    sim.add_component("plummer", np.array([mass, scale_radius]))


def test_clear_restores_defaults_and_deallocates_arrays():
    _set_basic_problem(n_particles=16, nsteps=20)
    sim.trim_orbits(2)
    sim.finalize()

    assert sim.x is not None
    assert sim.timestamps is not None
    assert sim.orbits is not None

    sim.max_ram_mb = 256.0
    sim.orbit_ram_limit_mb = 128.0
    sim.clear()

    assert sim.nparticles == 0
    assert sim.x is None
    assert sim.timestamps is None
    assert sim.orbits is None
    assert sim.max_ram_mb == pytest.approx(DEFAULT_MB)
    assert sim.orbit_ram_limit_mb == pytest.approx(DEFAULT_MB)


def test_setinitialconditions_memory_guard_blocks_oversized_setup():
    sim.max_ram_mb = 1.0e-9

    n_particles = 10
    arr = np.ones(n_particles)
    sim.setinitialconditions(arr, arr, arr, arr, arr, arr)

    assert sim.nparticles == 0
    assert sim.x is None


def test_finalize_builds_timestamps_and_allocates_orbits():
    n_particles = 24
    nsteps = 100
    nskip = 5

    _set_basic_problem(n_particles=n_particles, nsteps=nsteps)
    sim.trim_orbits(nskip)
    sim.finalize()

    assert sim.timestamps.shape == (nsteps + 1,)

    expected_saved_steps = nsteps // nskip + 1
    assert sim.orbits.shape == (expected_saved_steps, 6, n_particles)


def test_orbit_memory_guard_can_reduce_allocated_particles_to_zero():
    _set_basic_problem(n_particles=1000, nsteps=100)
    sim.orbit_ram_limit_mb = 1.0e-9
    sim.trim_orbits(1)
    sim.finalize()

    assert sim.orbits.shape[2] == 0


def test_finalize_uses_user_provided_timestamps_when_set():
    _set_basic_problem(n_particles=8, nsteps=200)
    user_timestamps = np.array([0.0, 0.5, 1.0, 1.5], dtype=float)

    sim.settimestamps(user_timestamps)
    sim.finalize()

    assert np.allclose(sim.timestamps, user_timestamps)




def test_run_without_gravity_refuses_to_start():
    n_particles = 8
    nsteps = 5
    dt = 0.1

    x0 = np.arange(n_particles, dtype=float)
    vx0 = np.ones(n_particles) * 2.0
    zeros = np.zeros(n_particles)

    sim.setinitialconditions(x0, zeros, zeros, vx0, zeros, zeros)
    sim.setscheme("leapfrog", [0.0, dt, nsteps])
    sim.run()

    assert not sim.run_success
    assert np.allclose(sim.x, x0)
    assert np.allclose(sim.vx, vx0)



def test_run_orbit_snapshot_count_matches_nskip():
    """Orbit array first dimension should equal nsteps//nskip + 1 (including t=0)."""
    n_particles = 5
    nsteps = 100
    nskip = 10

    _set_basic_problem(n_particles=n_particles, nsteps=nsteps)
    sim.trim_orbits(nskip)
    sim.run()

    expected_snapshots = nsteps // nskip + 1
    assert sim.orbits.shape[0] == expected_snapshots
