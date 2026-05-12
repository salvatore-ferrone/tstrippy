"""
Tests for the gravity module lifecycle API.
Covers: cleargravity, setgravityconstant, addgravitycomponent,
        finalizegravity, force, potential,
        printgravitystate.
"""
import numpy as np
import pytest
import tstrippy


@pytest.fixture(autouse=True)
def reset_gravity():
    """Always start each test with a clean gravity state."""
    tstrippy.gravity.gravity_clear()
    yield
    tstrippy.gravity.gravity_clear()


# ---------------------------------------------------------------------------
# Module state after clear
# ---------------------------------------------------------------------------

def test_cleargravity_resets_ncomp():
    g = tstrippy.gravity
    g.gravity_add_component("plummer", np.array([1e11, 2.0]))
    g.gravity_clear()
    assert g.gravity_ncomp == 0


def test_cleargravity_resets_finalized():
    g = tstrippy.gravity
    g.gravity_add_component("plummer", np.array([1e11, 2.0]))
    g.gravity_finalize()
    g.gravity_clear()
    assert not g.gravity_finalized


def test_cleargravity_resets_g_to_default():
    g = tstrippy.gravity
    g.gravity_set_gravity_constant(1.0)
    g.gravity_clear()
    assert g.gravity_g == pytest.approx(g.gravity_g_default, rel=1e-10)


# ---------------------------------------------------------------------------
# setgravityconstant
# ---------------------------------------------------------------------------

def test_setgravityconstant_stores_value():
    g = tstrippy.gravity
    g.gravity_set_gravity_constant(1.0)
    assert g.gravity_g == pytest.approx(1.0)


def test_setgravityconstant_marks_not_default():
    g = tstrippy.gravity
    g.gravity_set_gravity_constant(1.0)
    assert not g.gravity_g_is_default


def test_default_g_is_flagged():
    g = tstrippy.gravity
    assert g.gravity_g_is_default


# ---------------------------------------------------------------------------
# addgravitycomponent
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("model,params", [
    ("plummer",       np.array([1e11, 2.0])),
    ("hernquist",     np.array([5e10, 1.5])),
    ("miyamotonagai", np.array([6e10, 3.0, 0.3])),
    ("longmuralibar", np.array([1e10, 4.0, 0.4, 0.25])),
])
def test_addgravitycomponent_valid_models(model, params):
    g = tstrippy.gravity
    g.gravity_add_component(model, params)
    assert g.gravity_ncomp == 1


def test_addgravitycomponent_increments_ncomp():
    g = tstrippy.gravity
    g.gravity_add_component("plummer",   np.array([1e11, 2.0]))
    g.gravity_add_component("hernquist", np.array([5e10, 1.5]))
    assert g.gravity_ncomp == 2


def test_addgravitycomponent_invalid_model_stops():
    g = tstrippy.gravity
    g.gravity_add_component("bogusmodel", np.array([1.0, 2.0]))
    # Invalid model should NOT be added (ncomp unchanged)
    assert g.gravity_ncomp == 0


def test_addgravitycomponent_wrong_nparams_stops():
    g = tstrippy.gravity
    # plummer needs 2, we pass 3
    g.gravity_add_component("plummer", np.array([1e11, 2.0, 99.0]))
    # Wrong nparams should NOT be added (ncomp unchanged)
    assert g.gravity_ncomp == 0


def test_addgravitycomponent_after_finalize_stops():
    g = tstrippy.gravity
    g.gravity_add_component("plummer", np.array([1e11, 2.0]))
    g.gravity_finalize()
    ncomp_before = g.gravity_ncomp
    g.gravity_add_component("hernquist", np.array([5e10, 1.5]))
    # Adding after finalize should NOT modify ncomp
    assert g.gravity_ncomp == ncomp_before


# ---------------------------------------------------------------------------
# finalizegravity
# ---------------------------------------------------------------------------

def test_finalizegravity_sets_flag():
    g = tstrippy.gravity
    g.gravity_add_component("plummer", np.array([1e11, 2.0]))
    g.gravity_finalize()
    assert g.gravity_finalized


def test_finalizegravity_no_components_stops():
    g = tstrippy.gravity
    g.gravity_finalize()
    # Finalizing with no components should NOT set the finalized flag
    assert not g.gravity_finalized


# ---------------------------------------------------------------------------
# force
# ---------------------------------------------------------------------------

def test_force_returns_arrays():
    g = tstrippy.gravity
    g.gravity_add_component("plummer", np.array([1e11, 2.0]))
    g.gravity_finalize()
    x = np.array([8.0]);  y = np.zeros(1);  z = np.zeros(1)
    ax, ay, az = g.gravity_eval_force(x, y, z)
    assert ax.shape == (1,)
    assert ay.shape == (1,)
    assert az.shape == (1,)



def test_force_direction_on_axis():
    """Force on +x side of a spherical potential must point in -x."""
    g = tstrippy.gravity
    g.gravity_add_component("plummer", np.array([1e11, 2.0]))
    g.gravity_finalize()
    x = np.array([8.0]);  y = np.zeros(1);  z = np.zeros(1)
    ax, ay, az = g.gravity_eval_force(x, y, z)
    assert ax[0] < 0.0
    assert abs(ay[0]) < 1e-10
    assert abs(az[0]) < 1e-10


def test_force_superposition():
    """Two-component force is sum of individual forces."""
    g = tstrippy.gravity
    x = np.array([8.0]);  y = np.zeros(1);  z = np.zeros(1)

    g.gravity_add_component("plummer",   np.array([1e11, 2.0]))
    g.gravity_finalize()
    ax1, _, _ = g.gravity_eval_force(x, y, z)

    g.gravity_clear()
    g.gravity_add_component("hernquist", np.array([5e10, 1.5]))
    g.gravity_finalize()
    ax2, _, _ = g.gravity_eval_force(x, y, z)

    g.gravity_clear()
    g.gravity_add_component("plummer",   np.array([1e11, 2.0]))
    g.gravity_add_component("hernquist", np.array([5e10, 1.5]))
    g.gravity_finalize()
    ax_both, _, _ = g.gravity_eval_force(x, y, z)

    assert ax_both[0] == pytest.approx(ax1[0] + ax2[0], rel=1e-10)


# ---------------------------------------------------------------------------
# potential
# ---------------------------------------------------------------------------

def test_potential_returns_array():
    g = tstrippy.gravity
    g.gravity_add_component("plummer", np.array([1e11, 2.0]))
    g.gravity_finalize()
    x = np.array([8.0]);  y = np.zeros(1);  z = np.zeros(1)
    phi = g.gravity_eval_potential(x, y, z)
    assert phi.shape == (1,)


def test_potential_negative():
    """Gravitational potential is negative everywhere outside the source."""
    g = tstrippy.gravity
    g.gravity_add_component("plummer", np.array([1e11, 2.0]))
    g.gravity_finalize()
    x = np.array([8.0, 20.0]);  y = np.zeros(2);  z = np.zeros(2)
    phi = g.gravity_eval_potential(x, y, z)
    assert np.all(phi < 0.0)


def test_potential_decreases_with_distance():
    """Potential becomes less negative (increases) as distance grows."""
    g = tstrippy.gravity
    g.gravity_add_component("plummer", np.array([1e11, 2.0]))
    g.gravity_finalize()
    x = np.array([5.0, 10.0, 20.0]);  y = np.zeros(3);  z = np.zeros(3)
    phi = g.gravity_eval_potential(x, y, z)
    assert phi[0] < phi[1] < phi[2]



# ---------------------------------------------------------------------------
# Non-analytic lifecycle components
# ---------------------------------------------------------------------------

def test_exponential_oblate_halo_force_and_potential_are_finite():
    g = tstrippy.gravity
    g.gravity_add_component("exponential_oblate_halo", np.array([1.0e8, 5.0, 0.8]))
    g.gravity_finalize()

    x = np.array([8.0]); y = np.array([0.0]); z = np.array([0.5])
    ax, ay, az = g.gravity_eval_force(x, y, z)
    phi = g.gravity_eval_potential(x, y, z)

    assert np.isfinite(ax[0])
    assert np.isfinite(ay[0])
    assert np.isfinite(az[0])
    assert np.isfinite(phi[0])


def test_exponential_disk_bessel_phi_even_and_az_odd_in_z():
    g = tstrippy.gravity
    g.gravity_add_component("exponential_disk_bessel", np.array([5.0e8, 3.0, 0.3]))
    g.gravity_finalize()

    x = np.array([8.0, 8.0])
    y = np.zeros(2)
    z = np.array([0.4, -0.4])

    ax, ay, az = g.gravity_eval_force(x, y, z)
    phi = g.gravity_eval_potential(x, y, z)

    assert np.all(np.isfinite(ax))
    assert np.all(np.isfinite(ay))
    assert np.all(np.isfinite(az))
    assert np.all(np.isfinite(phi))
    assert phi[0] == pytest.approx(phi[1], rel=5e-4, abs=5e-7)
    assert az[0] == pytest.approx(-az[1], rel=5e-4, abs=5e-7)
