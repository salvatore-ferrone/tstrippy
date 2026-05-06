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


def test_finalize_precomputes_sh_for_single_component():
    g = gravitymini.gravity
    sh = gravitymini.sphericalharmonicsbfe

    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.finalizegravity()

    assert bool(sh.basis_expansion_initialized)


def test_finalize_initializes_component_phi_storage_for_multi_sh():
    g = gravitymini.gravity
    sh = gravitymini.sphericalharmonicsbfe

    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.addgravitycomponent("exponentialoblatehalo", [1.0, 1.0, 1.0])
    g.finalizegravity()

    assert int(g.gravity_ncomp) == 2
    assert int(sh.basis_ncomp) == 2
    assert sh.basis_phi_l_component_grid.shape[2] == 2


def test_store_component_phi_copies_active_phi_table():
    sh = gravitymini.sphericalharmonicsbfe

    sh.defaultinitsphericalharmonicbasis()
    sh.initsphericalharmoniccomponentphi(2)

    sh.basis_phi_l_grid[:, :] = 3.25
    sh.storesphericalharmoniccomponentphi(2)

    np.testing.assert_allclose(sh.basis_phi_l_component_grid[:, :, 1], sh.basis_phi_l_grid)


def test_finalize_populates_each_multi_sh_component_slot():
    g = gravitymini.gravity
    sh = gravitymini.sphericalharmonicsbfe

    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.addgravitycomponent("exponentialoblatehalo", [1.0, 1.0, 1.0])
    g.finalizegravity()

    slot1 = sh.basis_phi_l_component_grid[:, :, 0]
    slot2 = sh.basis_phi_l_component_grid[:, :, 1]

    assert np.any(np.abs(slot1) > 0.0)
    assert np.any(np.abs(slot2) > 0.0)


def test_load_component_phi_restores_slot_into_active_phi_grid():
    sh = gravitymini.sphericalharmonicsbfe

    sh.defaultinitsphericalharmonicbasis()
    sh.initsphericalharmoniccomponentphi(2)

    sh.basis_phi_l_grid[:, :] = 1.25
    sh.storesphericalharmoniccomponentphi(1)
    sh.basis_phi_l_grid[:, :] = 2.50
    sh.storesphericalharmoniccomponentphi(2)

    sh.basis_phi_l_grid[:, :] = -9.0
    sh.loadsphericalharmoniccomponentphi(2)

    np.testing.assert_allclose(sh.basis_phi_l_grid, sh.basis_phi_l_component_grid[:, :, 1])


def test_multi_sh_potential_reads_from_stored_component_slots():
    g = gravitymini.gravity
    sh = gravitymini.sphericalharmonicsbfe

    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.addgravitycomponent("exponentialoblatehalo", [1.0, 1.0, 1.0])
    g.finalizegravity()

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.1, 0.0], dtype=float)
    z = np.array([0.2, 0.0, -0.3], dtype=float)

    phi_before = g.evaluategravitypotential(x, y, z)
    sh.basis_phi_l_component_grid[:, :, 0] = 0.0
    phi_after = g.evaluategravitypotential(x, y, z)

    assert not np.allclose(phi_before, phi_after, rtol=1e-12, atol=1e-12)


def test_multi_sh_force_components_read_from_stored_component_slots():
    g = gravitymini.gravity
    sh = gravitymini.sphericalharmonicsbfe

    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.addgravitycomponent("exponentialoblatehalo", [1.0, 1.0, 1.0])
    g.finalizegravity()

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.1, 0.0], dtype=float)
    z = np.array([0.2, 0.0, -0.3], dtype=float)

    ax_before, ay_before, az_before = g.evaluategravityforcecomponents(x, y, z)
    sh.basis_phi_l_component_grid[:, :, 0] = 0.0
    ax_after, ay_after, az_after = g.evaluategravityforcecomponents(x, y, z)

    before = np.concatenate((ax_before, ay_before, az_before), axis=0)
    after = np.concatenate((ax_after, ay_after, az_after), axis=0)
    assert not np.allclose(before, after, rtol=1e-12, atol=1e-12)


def test_multi_sh_total_force_reads_from_stored_component_slots():
    g = gravitymini.gravity
    sh = gravitymini.sphericalharmonicsbfe

    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.addgravitycomponent("exponentialoblatehalo", [1.0, 1.0, 1.0])
    g.finalizegravity()

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.1, 0.0], dtype=float)
    z = np.array([0.2, 0.0, -0.3], dtype=float)

    ax_before, ay_before, az_before = g.evaluategravityforces(x, y, z)
    sh.basis_phi_l_component_grid[:, :, 0] = 0.0
    ax_after, ay_after, az_after = g.evaluategravityforces(x, y, z)

    before = np.concatenate((ax_before, ay_before, az_before), axis=0)
    after = np.concatenate((ax_after, ay_after, az_after), axis=0)
    assert not np.allclose(before, after, rtol=1e-12, atol=1e-12)


def test_interleaved_components_preserve_sh_slot_mapping_order():
    g = gravitymini.gravity
    sh = gravitymini.sphericalharmonicsbfe

    g.cleargravity()
    g.addgravitycomponent("plummer", [1.0, 1.0])
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.addgravitycomponent("hernquist", [1.0, 1.0])
    g.addgravitycomponent("exponentialoblatehalo", [1.0, 1.0, 1.0])
    g.finalizegravity()

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.1, 0.0], dtype=float)
    z = np.array([0.2, 0.0, -0.3], dtype=float)

    phi0 = g.evaluategravitypotential(x, y, z)

    slot1_backup = sh.basis_phi_l_component_grid[:, :, 0].copy()
    sh.basis_phi_l_component_grid[:, :, 0] = 0.0
    phi1 = g.evaluategravitypotential(x, y, z)
    sh.basis_phi_l_component_grid[:, :, 0] = slot1_backup

    assert not np.allclose(phi0, phi1, rtol=1e-12, atol=1e-12)

    slot2_backup = sh.basis_phi_l_component_grid[:, :, 1].copy()
    sh.basis_phi_l_component_grid[:, :, 1] = 0.0
    phi2 = g.evaluategravitypotential(x, y, z)
    sh.basis_phi_l_component_grid[:, :, 1] = slot2_backup

    assert not np.allclose(phi0, phi2, rtol=1e-12, atol=1e-12)
    assert not np.allclose(phi1, phi2, rtol=1e-12, atol=1e-12)