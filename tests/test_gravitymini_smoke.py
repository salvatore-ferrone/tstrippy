import numpy as np
import tstrippy


def test_ibata_density_is_exposed():
    g = tstrippy.gravity
    assert hasattr(g, "ibata2024halo_density")


def _reset_with_single_plummer():
    g = tstrippy.gravity
    g.cleargravity()
    g.addgravitycomponent("plummer", [1.0, 1.0])
    g.finalizegravity()
    return g


def test_component_api_returns_fixed_shape():
    g = _reset_with_single_plummer()

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.0, 0.0], dtype=float)
    z = np.array([0.0, 0.0, 0.0], dtype=float)

    ax_comp, ay_comp, az_comp = g.force_components(x, y, z)

    assert ax_comp.shape == (16, 3)
    assert ay_comp.shape == (16, 3)
    assert az_comp.shape == (16, 3)


def test_component_sum_matches_total_force():
    g = _reset_with_single_plummer()

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.0, 0.0], dtype=float)
    z = np.array([0.0, 0.0, 0.0], dtype=float)

    ax_comp, ay_comp, az_comp = g.force_components(x, y, z)
    ax, ay, az = g.force(x, y, z)

    ncomp = int(g.gravity_ncomp)
    np.testing.assert_allclose(ax_comp[:ncomp, :].sum(axis=0), ax, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(ay_comp[:ncomp, :].sum(axis=0), ay, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(az_comp[:ncomp, :].sum(axis=0), az, rtol=1e-12, atol=1e-12)


def test_ibata2024halo():
    g = tstrippy.gravity

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.0, 0.0], dtype=float)
    z = np.array([0.0, 0.0, 0.0], dtype=float)


    params = [1, 1, 100, .8, 1.4, 3.0]
    den = g.ibata2024halo_density(params,x,y,z)
    
    assert np.all(den) > 0 


def test_ibata_component_force_api_runs():
    g = tstrippy.gravity
    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.finalizegravity()

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.1, 0.0], dtype=float)
    z = np.array([0.2, 0.0, -0.3], dtype=float)

    ax_comp, ay_comp, az_comp = g.force_components(x, y, z)

    assert np.all(np.isfinite(ax_comp[0, :]))
    assert np.all(np.isfinite(ay_comp[0, :]))
    assert np.all(np.isfinite(az_comp[0, :]))
    assert np.any(np.abs(ax_comp[0, :]) > 0.0) or np.any(np.abs(ay_comp[0, :]) > 0.0) or np.any(np.abs(az_comp[0, :]) > 0.0)


def test_ibata_total_potential_api_runs():
    g = tstrippy.gravity
    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.finalizegravity()

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.1, 0.0], dtype=float)
    z = np.array([0.2, 0.0, -0.3], dtype=float)

    phi = g.potential(x, y, z)

    assert np.all(np.isfinite(phi))
    assert np.any(np.abs(phi) > 0.0)


def test_finalize_precomputes_sh_for_single_component():
    g = tstrippy.gravity
    sh = tstrippy.sphericalharmonicsbfe

    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.finalizegravity()

    assert bool(sh.basis_expansion_initialized)


def test_finalize_initializes_component_phi_storage_for_multi_sh():
    g = tstrippy.gravity
    sh = tstrippy.sphericalharmonicsbfe

    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.addgravitycomponent("exponentialoblatehalo", [1.0, 1.0, 1.0])
    g.finalizegravity()

    assert int(g.gravity_ncomp) == 2
    assert int(sh.basis_ncomp) == 2
    assert sh.basis_phi_l_component_grid.shape[2] == 2


def test_store_component_phi_copies_active_phi_table():
    sh = tstrippy.sphericalharmonicsbfe

    sh.defaultinitsphericalharmonicbasis()
    sh.initsphericalharmoniccomponentphi(2)

    sh.basis_phi_l_grid[:, :] = 3.25
    sh.storesphericalharmoniccomponentphi(2)

    np.testing.assert_allclose(sh.basis_phi_l_component_grid[:, :, 1], sh.basis_phi_l_grid)


def test_finalize_populates_each_multi_sh_component_slot():
    g = tstrippy.gravity
    sh = tstrippy.sphericalharmonicsbfe

    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.addgravitycomponent("exponentialoblatehalo", [1.0, 1.0, 1.0])
    g.finalizegravity()

    slot1 = sh.basis_phi_l_component_grid[:, :, 0]
    slot2 = sh.basis_phi_l_component_grid[:, :, 1]

    assert np.any(np.abs(slot1) > 0.0)
    assert np.any(np.abs(slot2) > 0.0)


def test_load_component_phi_restores_slot_into_active_phi_grid():
    sh = tstrippy.sphericalharmonicsbfe

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
    g = tstrippy.gravity
    sh = tstrippy.sphericalharmonicsbfe

    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.addgravitycomponent("exponentialoblatehalo", [1.0, 1.0, 1.0])
    g.finalizegravity()

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.1, 0.0], dtype=float)
    z = np.array([0.2, 0.0, -0.3], dtype=float)

    phi_before = g.potential(x, y, z)
    sh.basis_phi_l_component_grid[:, :, 0] = 0.0
    phi_after = g.potential(x, y, z)

    assert not np.allclose(phi_before, phi_after, rtol=1e-12, atol=1e-12)


def test_multi_sh_force_components_read_from_stored_component_slots():
    g = tstrippy.gravity
    sh = tstrippy.sphericalharmonicsbfe

    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.addgravitycomponent("exponentialoblatehalo", [1.0, 1.0, 1.0])
    g.finalizegravity()

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.1, 0.0], dtype=float)
    z = np.array([0.2, 0.0, -0.3], dtype=float)

    ax_before, ay_before, az_before = g.force_components(x, y, z)
    sh.basis_phi_l_component_grid[:, :, 0] = 0.0
    ax_after, ay_after, az_after = g.force_components(x, y, z)

    before = np.concatenate((ax_before, ay_before, az_before), axis=0)
    after = np.concatenate((ax_after, ay_after, az_after), axis=0)
    assert not np.allclose(before, after, rtol=1e-12, atol=1e-12)


def test_multi_sh_total_force_reads_from_stored_component_slots():
    g = tstrippy.gravity
    sh = tstrippy.sphericalharmonicsbfe

    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.addgravitycomponent("exponentialoblatehalo", [1.0, 1.0, 1.0])
    g.finalizegravity()

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.1, 0.0], dtype=float)
    z = np.array([0.2, 0.0, -0.3], dtype=float)

    ax_before, ay_before, az_before = g.force(x, y, z)
    sh.basis_phi_l_component_grid[:, :, 0] = 0.0
    ax_after, ay_after, az_after = g.force(x, y, z)

    before = np.concatenate((ax_before, ay_before, az_before), axis=0)
    after = np.concatenate((ax_after, ay_after, az_after), axis=0)
    assert not np.allclose(before, after, rtol=1e-12, atol=1e-12)


def test_interleaved_components_preserve_sh_slot_mapping_order():
    g = tstrippy.gravity
    sh = tstrippy.sphericalharmonicsbfe

    g.cleargravity()
    g.addgravitycomponent("plummer", [1.0, 1.0])
    g.addgravitycomponent("ibata2024halo", [1.0, 1.0, 100.0, 0.8, 1.4, 3.0])
    g.addgravitycomponent("hernquist", [1.0, 1.0])
    g.addgravitycomponent("exponentialoblatehalo", [1.0, 1.0, 1.0])
    g.finalizegravity()

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.1, 0.0], dtype=float)
    z = np.array([0.2, 0.0, -0.3], dtype=float)

    phi0 = g.potential(x, y, z)

    slot1_backup = sh.basis_phi_l_component_grid[:, :, 0].copy()
    sh.basis_phi_l_component_grid[:, :, 0] = 0.0
    phi1 = g.potential(x, y, z)
    sh.basis_phi_l_component_grid[:, :, 0] = slot1_backup

    assert not np.allclose(phi0, phi1, rtol=1e-12, atol=1e-12)

    slot2_backup = sh.basis_phi_l_component_grid[:, :, 1].copy()
    sh.basis_phi_l_component_grid[:, :, 1] = 0.0
    phi2 = g.potential(x, y, z)
    sh.basis_phi_l_component_grid[:, :, 1] = slot2_backup

    assert not np.allclose(phi0, phi2, rtol=1e-12, atol=1e-12)
    assert not np.allclose(phi1, phi2, rtol=1e-12, atol=1e-12)


def test_ibata_exponential_force_components_match_independent_components():
    g = tstrippy.gravity

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.1, 0.0], dtype=float)
    z = np.array([0.2, 0.0, -0.3], dtype=float)

    ibata = [1.0, 1.0, 100.0, 0.8, 1.4, 3.0]
    expo = [1.0, 1.0, 1.0]

    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", ibata)
    g.finalizegravity()
    ax_i, ay_i, az_i = g.force(x, y, z)

    g.cleargravity()
    g.addgravitycomponent("exponentialoblatehalo", expo)
    g.finalizegravity()
    ax_e, ay_e, az_e = g.force(x, y, z)

    g.cleargravity()
    g.addgravitycomponent("ibata2024halo", ibata)
    g.addgravitycomponent("exponentialoblatehalo", expo)
    g.finalizegravity()
    ax_c, ay_c, az_c = g.force_components(x, y, z)

    np.testing.assert_allclose(ax_c[0, :], ax_i, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(ay_c[0, :], ay_i, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(az_c[0, :], az_i, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(ax_c[1, :], ax_e, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(ay_c[1, :], ay_e, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(az_c[1, :], az_e, rtol=1e-12, atol=1e-12)


def test_all_component_pairs_are_commutative_and_force_consistent_on_small_grid():
    g = tstrippy.gravity

    def rms(a, b):
        return np.sqrt(np.mean((a.ravel() - b.ravel()) ** 2))

    rho_h_table = 11.4
    rho0_halo = (rho_h_table / 1000.0) * 1e9
    r0 = 14.7
    rt = 1e3
    q = 0.5
    gamma = 1.0
    beta = 3.0

    components = [
        ("plummer", [1e12, r0]),
        ("ibata2024halo", [rho0_halo, r0, rt, q, gamma, beta]),
        ("exponentialoblatehalo", [rho0_halo / 2.0, r0 / 2.0, q]),
        ("hernquist", [1e10, r0]),
    ]

    # Compact diagnostics grid: 10x10 points over +/-2*scale radius in x-z plane.
    x = np.linspace(-2.0 * r0, 2.0 * r0, 10)
    z = np.linspace(-2.0 * r0, 2.0 * r0, 10)
    X, Z = np.meshgrid(x, z, indexing="xy")
    xf = X.ravel()
    zf = Z.ravel()
    y0 = np.zeros_like(xf)

    single_force = {}
    for name, params in components:
        g.cleargravity()
        g.addgravitycomponent(name, params)
        g.finalizegravity()
        single_force[name] = g.force(xf, y0, zf)

    for name_i, params_i in components:
        for name_j, params_j in components:
            g.cleargravity()
            g.addgravitycomponent(name_i, params_i)
            g.addgravitycomponent(name_j, params_j)
            g.finalizegravity()

            ax, ay, az = g.force(xf, y0, zf)
            phi = g.potential(xf, y0, zf)
            ax_c, ay_c, az_c = g.force_components(xf, y0, zf)

            ax_sum = ax_c.sum(axis=0)
            ay_sum = ay_c.sum(axis=0)
            az_sum = az_c.sum(axis=0)

            g.cleargravity()
            g.addgravitycomponent(name_j, params_j)
            g.addgravitycomponent(name_i, params_i)
            g.finalizegravity()
            ax_r, ay_r, az_r = g.force(xf, y0, zf)
            phi_r = g.potential(xf, y0, zf)

            ax_s = single_force[name_i][0] + single_force[name_j][0]
            ay_s = single_force[name_i][1] + single_force[name_j][1]
            az_s = single_force[name_i][2] + single_force[name_j][2]

            assert rms(ax, ax_sum) < 1e-12
            assert rms(ay, ay_sum) < 1e-12
            assert rms(az, az_sum) < 1e-12

            assert rms(ax, ax_s) < 1e-12
            assert rms(ay, ay_s) < 1e-12
            assert rms(az, az_s) < 1e-12

            assert rms(ax, ax_r) < 1e-12
            assert rms(ay, ay_r) < 1e-12
            assert rms(az, az_r) < 1e-12
            assert rms(phi, phi_r) < 1e-12


def test_component_potential_sum_matches_total_potential():
    g = tstrippy.gravity

    g.cleargravity()
    g.addgravitycomponent("plummer", [1e12, 14.7])
    g.addgravitycomponent("ibata2024halo", [1.14e7, 14.7, 1e3, 0.5, 1.0, 3.0])
    g.addgravitycomponent("exponentialoblatehalo", [5.7e6, 7.35, 0.5])
    g.finalizegravity()

    x = np.array([1.0, 2.0, 3.0], dtype=float)
    y = np.array([0.0, 0.1, 0.0], dtype=float)
    z = np.array([0.2, 0.0, -0.3], dtype=float)

    phi_total = g.potential(x, y, z)
    phi_comp = g.potential_components(x, y, z)

    ncomp = int(g.gravity_ncomp)
    np.testing.assert_allclose(phi_comp[:ncomp, :].sum(axis=0), phi_total, rtol=1e-12, atol=1e-12)