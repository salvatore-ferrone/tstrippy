"""
Tests for the gravity module lifecycle API.
Covers: cleargravity, setgravityconstant, addgravitycomponent,
        finalizegravity, evaluategravityforces, evaluategravitypotential,
        printgravitystate.
"""
import numpy as np
import pytest
import tstrippy


@pytest.fixture(autouse=True)
def reset_gravity():
    """Always start each test with a clean gravity state."""
    tstrippy.gravity.cleargravity()
    yield
    tstrippy.gravity.cleargravity()


# ---------------------------------------------------------------------------
# Module state after clear
# ---------------------------------------------------------------------------

def test_cleargravity_resets_ncomp():
    g = tstrippy.gravity
    g.addgravitycomponent("plummer", np.array([1e11, 2.0]))
    g.cleargravity()
    assert g.gravity_ncomp == 0


def test_cleargravity_resets_finalized():
    g = tstrippy.gravity
    g.addgravitycomponent("plummer", np.array([1e11, 2.0]))
    g.finalizegravity()
    g.cleargravity()
    assert not g.gravity_finalized


def test_cleargravity_resets_g_to_default():
    g = tstrippy.gravity
    g.setgravityconstant(1.0)
    g.cleargravity()
    assert g.gravity_g == pytest.approx(g.gravity_g_default, rel=1e-10)


# ---------------------------------------------------------------------------
# setgravityconstant
# ---------------------------------------------------------------------------

def test_setgravityconstant_stores_value():
    g = tstrippy.gravity
    g.setgravityconstant(1.0)
    assert g.gravity_g == pytest.approx(1.0)


def test_setgravityconstant_marks_not_default():
    g = tstrippy.gravity
    g.setgravityconstant(1.0)
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
    g.addgravitycomponent(model, params)
    assert g.gravity_ncomp == 1


def test_addgravitycomponent_increments_ncomp():
    g = tstrippy.gravity
    g.addgravitycomponent("plummer",   np.array([1e11, 2.0]))
    g.addgravitycomponent("hernquist", np.array([5e10, 1.5]))
    assert g.gravity_ncomp == 2


def test_addgravitycomponent_invalid_model_stops():
    g = tstrippy.gravity
    g.addgravitycomponent("bogusmodel", np.array([1.0, 2.0]))
    # Invalid model should NOT be added (ncomp unchanged)
    assert g.gravity_ncomp == 0


def test_addgravitycomponent_wrong_nparams_stops():
    g = tstrippy.gravity
    # plummer needs 2, we pass 3
    g.addgravitycomponent("plummer", np.array([1e11, 2.0, 99.0]))
    # Wrong nparams should NOT be added (ncomp unchanged)
    assert g.gravity_ncomp == 0


def test_addgravitycomponent_after_finalize_stops():
    g = tstrippy.gravity
    g.addgravitycomponent("plummer", np.array([1e11, 2.0]))
    g.finalizegravity()
    ncomp_before = g.gravity_ncomp
    g.addgravitycomponent("hernquist", np.array([5e10, 1.5]))
    # Adding after finalize should NOT modify ncomp
    assert g.gravity_ncomp == ncomp_before


# ---------------------------------------------------------------------------
# finalizegravity
# ---------------------------------------------------------------------------

def test_finalizegravity_sets_flag():
    g = tstrippy.gravity
    g.addgravitycomponent("plummer", np.array([1e11, 2.0]))
    g.finalizegravity()
    assert g.gravity_finalized


def test_finalizegravity_no_components_stops():
    g = tstrippy.gravity
    g.finalizegravity()
    # Finalizing with no components should NOT set the finalized flag
    assert not g.gravity_finalized


# ---------------------------------------------------------------------------
# evaluategravityforces
# ---------------------------------------------------------------------------

def test_evaluategravityforces_returns_arrays():
    g = tstrippy.gravity
    g.addgravitycomponent("plummer", np.array([1e11, 2.0]))
    g.finalizegravity()
    x = np.array([8.0]);  y = np.zeros(1);  z = np.zeros(1)
    ax, ay, az = g.evaluategravityforces(x, y, z)
    assert ax.shape == (1,)
    assert ay.shape == (1,)
    assert az.shape == (1,)


def test_evaluategravityforces_before_finalize_autofinalizes_for_analytic():
    g = tstrippy.gravity
    g.addgravitycomponent("plummer", np.array([1e11, 2.0]))
    x = np.array([8.0]);  y = np.zeros(1);  z = np.zeros(1)
    ax, ay, az = g.evaluategravityforces(x, y, z)
    # Analytic-only configurations auto-finalize on first evaluate call.
    assert g.gravity_finalized
    assert ax[0] != 0.0
    assert abs(ay[0]) < 1e-10
    assert abs(az[0]) < 1e-10


def test_evaluategravityforces_direction_on_axis():
    """Force on +x side of a spherical potential must point in -x."""
    g = tstrippy.gravity
    g.addgravitycomponent("plummer", np.array([1e11, 2.0]))
    g.finalizegravity()
    x = np.array([8.0]);  y = np.zeros(1);  z = np.zeros(1)
    ax, ay, az = g.evaluategravityforces(x, y, z)
    assert ax[0] < 0.0
    assert abs(ay[0]) < 1e-10
    assert abs(az[0]) < 1e-10


def test_evaluategravityforces_superposition():
    """Two-component force is sum of individual forces."""
    g = tstrippy.gravity
    x = np.array([8.0]);  y = np.zeros(1);  z = np.zeros(1)

    g.addgravitycomponent("plummer",   np.array([1e11, 2.0]))
    g.finalizegravity()
    ax1, _, _ = g.evaluategravityforces(x, y, z)

    g.cleargravity()
    g.addgravitycomponent("hernquist", np.array([5e10, 1.5]))
    g.finalizegravity()
    ax2, _, _ = g.evaluategravityforces(x, y, z)

    g.cleargravity()
    g.addgravitycomponent("plummer",   np.array([1e11, 2.0]))
    g.addgravitycomponent("hernquist", np.array([5e10, 1.5]))
    g.finalizegravity()
    ax_both, _, _ = g.evaluategravityforces(x, y, z)

    assert ax_both[0] == pytest.approx(ax1[0] + ax2[0], rel=1e-10)


# ---------------------------------------------------------------------------
# evaluategravitypotential
# ---------------------------------------------------------------------------

def test_evaluategravitypotential_returns_array():
    g = tstrippy.gravity
    g.addgravitycomponent("plummer", np.array([1e11, 2.0]))
    g.finalizegravity()
    x = np.array([8.0]);  y = np.zeros(1);  z = np.zeros(1)
    phi = g.evaluategravitypotential(x, y, z)
    assert phi.shape == (1,)


def test_evaluategravitypotential_negative():
    """Gravitational potential is negative everywhere outside the source."""
    g = tstrippy.gravity
    g.addgravitycomponent("plummer", np.array([1e11, 2.0]))
    g.finalizegravity()
    x = np.array([8.0, 20.0]);  y = np.zeros(2);  z = np.zeros(2)
    phi = g.evaluategravitypotential(x, y, z)
    assert np.all(phi < 0.0)


def test_evaluategravitypotential_decreases_with_distance():
    """Potential becomes less negative (increases) as distance grows."""
    g = tstrippy.gravity
    g.addgravitycomponent("plummer", np.array([1e11, 2.0]))
    g.finalizegravity()
    x = np.array([5.0, 10.0, 20.0]);  y = np.zeros(3);  z = np.zeros(3)
    phi = g.evaluategravitypotential(x, y, z)
    assert phi[0] < phi[1] < phi[2]


def test_evaluategravitypotential_before_finalize_autofinalizes_for_analytic():
    g = tstrippy.gravity
    g.addgravitycomponent("plummer", np.array([1e11, 2.0]))
    x = np.array([8.0]);  y = np.zeros(1);  z = np.zeros(1)
    phi = g.evaluategravitypotential(x, y, z)
    # Analytic-only configurations auto-finalize on first evaluate call.
    assert g.gravity_finalized
    assert phi[0] < 0.0
