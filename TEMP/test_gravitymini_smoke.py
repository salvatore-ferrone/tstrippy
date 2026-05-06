import numpy as np
import gravitymini


def test_ibata_density_is_exposed():
    g = gravitymini.gravity
    assert hasattr(g, "ibata2024halo_density")


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


def test_ibata2024halo():
    g = gravitymini.gravity

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.0, 0.0], dtype=float)
    z = np.array([0.0, 0.0, 0.0], dtype=float)


    params = [1, 1, 100, .8, 1.4, 3.0]
    den = g.ibata2024halo_density(params,x,y,z)
    
    assert np.all(den) > 0 


def test_ibata_component_force_api_runs():
    g = gravitymini.gravity
    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.finalizegravity()

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.1, 0.0], dtype=float)
    z = np.array([0.2, 0.0, -0.3], dtype=float)

    ax_comp, ay_comp, az_comp = g.evaluategravityforcecomponents(x, y, z)

    assert np.all(np.isfinite(ax_comp[0, :]))
    assert np.all(np.isfinite(ay_comp[0, :]))
    assert np.all(np.isfinite(az_comp[0, :]))
    assert np.any(np.abs(ax_comp[0, :]) > 0.0) or np.any(np.abs(ay_comp[0, :]) > 0.0) or np.any(np.abs(az_comp[0, :]) > 0.0)


def test_ibata_total_potential_api_runs():
    g = gravitymini.gravity
    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.finalizegravity()

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.1, 0.0], dtype=float)
    z = np.array([0.2, 0.0, -0.3], dtype=float)

    phi = g.evaluategravitypotential(x, y, z)

    assert np.all(np.isfinite(phi))
    assert np.any(np.abs(phi) > 0.0)