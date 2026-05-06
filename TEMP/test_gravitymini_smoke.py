import numpy as np
import gravitymini


def _reset_with_single_plummer():
    g = gravitymini.gravity
    g.cleargravity()
    g.addgravitycomponent("plummer", [1.0, 1.0])
    g.finalizegravity()
    return g


def test_component_api_returns_fixed_shape():
    g = _reset_with_single_plummer()

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.0, 0.0], dtype=float)
    z = np.array([0.0, 0.0, 0.0], dtype=float)

    ax_comp, ay_comp, az_comp = g.evaluategravityforcecomponents(x, y, z)

    assert ax_comp.shape == (16, 3)
    assert ay_comp.shape == (16, 3)
    assert az_comp.shape == (16, 3)


def test_component_sum_matches_total_force():
    g = _reset_with_single_plummer()

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.0, 0.0], dtype=float)
    z = np.array([0.0, 0.0, 0.0], dtype=float)

    ax_comp, ay_comp, az_comp = g.evaluategravityforcecomponents(x, y, z)
    ax, ay, az = g.evaluategravityforces(x, y, z)

    ncomp = int(g.gravity_ncomp)
    np.testing.assert_allclose(ax_comp[:ncomp, :].sum(axis=0), ax, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(ay_comp[:ncomp, :].sum(axis=0), ay, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(az_comp[:ncomp, :].sum(axis=0), az, rtol=1e-12, atol=1e-12)
