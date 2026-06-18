import math

import numpy as np
import pytest

import tstrippy


def _configure_single_bessel_disk(sigma0=1.0, hR=4.0, hZ=0.8, nr=128, nz=128, nk=256):
    g = tstrippy.gravity

    g.clear()
    g.bessel_initialize(nr, nz, nk)
    g.add_component("exponentialdisk", [sigma0, hR, hZ])
    g.bessel_set_component_scales(1, hR, hZ)
    g.finalize()
    return g


def _disk_total_mass(sigma0, hR):
    return 2.0 * math.pi * sigma0 * hR**2


def _expdisk_density(R, z, sigma0, hR, hZ):
    return (sigma0 / (2.0 * hZ)) * np.exp(-R / hR) * np.exp(-np.abs(z) / hZ)


def _far_field_inferred_masses(g, radii, theta):
    gconst = float(tstrippy.gravity.gravity_g_default)

    x = radii * np.cos(theta)
    y = np.zeros_like(x)
    z = radii * np.sin(theta)

    phi = g.potential(x, y, z)
    ax, ay, az = g.force(x, y, z)

    pos = np.column_stack((x, y, z))
    rhat = pos / np.linalg.norm(pos, axis=1, keepdims=True)
    acc = np.column_stack((ax, ay, az))
    a_r = np.sum(acc * rhat, axis=1)

    m_phi = -radii * phi / gconst
    m_a = -(radii**2) * a_r / gconst
    return m_phi, m_a


def _poisson_residual_map(g, sigma0, hR, hZ):
    gconst = float(tstrippy.gravity.gravity_g_default)

    R = np.linspace(0.6, 8.0, 33)
    Z = np.linspace(-3.0, 3.0, 49)
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")

    x = RR.ravel()
    y = np.zeros_like(x)
    z = ZZ.ravel()

    phi = g.potential(x, y, z).reshape(RR.shape)

    dR = R[1] - R[0]
    dZ = Z[1] - Z[0]

    phi_rr = (phi[2:, 1:-1] - 2.0 * phi[1:-1, 1:-1] + phi[:-2, 1:-1]) / (dR**2)
    phi_r = (phi[2:, 1:-1] - phi[:-2, 1:-1]) / (2.0 * dR)
    phi_zz = (phi[1:-1, 2:] - 2.0 * phi[1:-1, 1:-1] + phi[1:-1, :-2]) / (dZ**2)

    R_inner = RR[1:-1, 1:-1]
    Z_inner = ZZ[1:-1, 1:-1]

    laplacian = phi_rr + (phi_r / R_inner) + phi_zz
    rho = _expdisk_density(R_inner, Z_inner, sigma0, hR, hZ)
    source = 4.0 * math.pi * gconst * rho

    residual = laplacian - source
    rel = np.abs(residual) / np.maximum(np.abs(source), 1.0e-12)

    return {
        "residual": residual,
        "rel": rel,
        "source": source,
    }


def _configure_two_bessel_disks(
    sigma0_1=2.0,
    hR_1=4.0,
    hZ_1=2.0,
    sigma0_2=1.0,
    hR_2=6.0,
    hZ_2=0.5,
    nr=256,
    nz=256,
    nk=256,
    scale_divisor=2.0,
):
    g = tstrippy.gravity

    g.clear()
    g.bessel_initialize(nr, nz, nk)
    g.add_component("exponentialdisk", [sigma0_1, hR_1, hZ_1])
    g.add_component("exponentialdisk", [sigma0_2, hR_2, hZ_2])
    g.bessel_set_component_scales(1, hR_1 / scale_divisor, hZ_1 / scale_divisor)
    g.bessel_set_component_scales(2, hR_2 / scale_divisor, hZ_2 / scale_divisor)
    g.finalize()
    return g


def _direct_reference_potential(R_eval, z_eval, sigma0, hR, hZ, nR=140, nZ=120, nphi=180):
    gconst = float(tstrippy.gravity.gravity_g_default)

    r_max = 12.0 * hR
    z_max = 10.0 * hZ

    r_edges = np.linspace(0.0, r_max, nR + 1)
    z_edges = np.linspace(-z_max, z_max, nZ + 1)
    r_mid = 0.5 * (r_edges[:-1] + r_edges[1:])
    z_mid = 0.5 * (z_edges[:-1] + z_edges[1:])
    dR = r_edges[1] - r_edges[0]
    dZ = z_edges[1] - z_edges[0]

    phi_grid = np.linspace(0.0, 2.0 * math.pi, nphi, endpoint=False)
    cosphi = np.cos(phi_grid)
    dphi = 2.0 * math.pi / nphi

    phi_total = 0.0
    for r_src in r_mid:
        density_r = (sigma0 / (2.0 * hZ)) * math.exp(-r_src / hR)
        for z_src in z_mid:
            rho = density_r * math.exp(-abs(z_src) / hZ)
            sep2 = R_eval**2 + r_src**2 - 2.0 * R_eval * r_src * cosphi + (z_eval - z_src) ** 2
            sep = np.sqrt(np.maximum(sep2, 1.0e-14))
            phi_ring = np.sum(1.0 / sep) * dphi
            phi_total += -gconst * rho * r_src * dR * dZ * phi_ring

    return phi_total


def test_bessel_disk_even_symmetry_and_midplane_parity():
    g = _configure_single_bessel_disk()

    x = np.array([2.0, 2.0, 2.0], dtype=float)
    y = np.array([0.0, 0.0, 0.0], dtype=float)
    z = np.array([0.6, -0.6, 0.0], dtype=float)

    ax, ay, az = g.force(x, y, z)
    phi = g.potential(x, y, z)

    np.testing.assert_allclose(phi[0], phi[1], rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(ax[0], ax[1], rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(ay[0], ay[1], rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(az[0], -az[1], rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(az[2], 0.0, atol=1e-10)


def test_bessel_force_matches_finite_difference_of_potential():
    g = _configure_single_bessel_disk()
    delta = 1.0e-3

    base_x = np.array([2.0], dtype=float)
    base_y = np.array([0.0], dtype=float)
    base_z = np.array([0.4], dtype=float)

    phi_x_plus = g.potential(base_x + delta, base_y, base_z)[0]
    phi_x_minus = g.potential(base_x - delta, base_y, base_z)[0]
    phi_z_plus = g.potential(base_x, base_y, base_z + delta)[0]
    phi_z_minus = g.potential(base_x, base_y, base_z - delta)[0]

    ax, ay, az = g.force(base_x, base_y, base_z)

    dphi_dx = (phi_x_plus - phi_x_minus) / (2.0 * delta)
    dphi_dz = (phi_z_plus - phi_z_minus) / (2.0 * delta)

    np.testing.assert_allclose(ax[0], -dphi_dx, rtol=5e-3, atol=5e-5)
    np.testing.assert_allclose(az[0], -dphi_dz, rtol=5e-3, atol=5e-5)
    np.testing.assert_allclose(ay[0], 0.0, atol=1e-12)


def test_bessel_far_field_matches_monopole_mass():
    sigma0 = 1.0
    hR = 4.0
    hZ = 0.8
    g = _configure_single_bessel_disk(sigma0=sigma0, hR=hR, hZ=hZ)

    x = np.array([80.0], dtype=float)
    y = np.array([0.0], dtype=float)
    z = np.array([0.0], dtype=float)

    phi = g.potential(x, y, z)[0]
    ax, ay, az = g.force(x, y, z)

    mass = _disk_total_mass(sigma0, hR)
    gconst = float(tstrippy.gravity.gravity_g_default)
    radius = x[0]
    phi_ref = -gconst * mass / radius
    ax_ref = -gconst * mass / radius**2

    np.testing.assert_allclose(phi, phi_ref, rtol=2.5e-1, atol=1e-6)
    np.testing.assert_allclose(ax[0], ax_ref, rtol=2.5e-1, atol=1e-6)
    np.testing.assert_allclose(ay[0], 0.0, atol=1e-12)
    np.testing.assert_allclose(az[0], 0.0, atol=1e-12)


def test_bessel_far_field_force_aligns_with_negative_rhat():
    sigma0 = 1.0
    hR = 4.0
    hZ = 0.8
    g = _configure_single_bessel_disk(sigma0=sigma0, hR=hR, hZ=hZ, nr=192, nz=160, nk=384)

    # Probe several increasing radii on a fixed oblique ray so both radial and
    # vertical force components are exercised.
    radii = np.array([20.0, 40.0, 80.0, 120.0], dtype=float)
    theta = 0.35  # polar angle from midplane in the R-z meridional plane
    x = radii * np.cos(theta)
    y = np.zeros_like(x)
    z = radii * np.sin(theta)

    ax, ay, az = g.force(x, y, z)
    force = np.column_stack((ax, ay, az))
    pos = np.column_stack((x, y, z))

    minus_rhat = -pos / np.linalg.norm(pos, axis=1, keepdims=True)
    force_norm = np.linalg.norm(force, axis=1)

    assert np.all(force_norm > 0.0)

    cos_align = np.sum(force * minus_rhat, axis=1) / force_norm

    # Must be attractive everywhere along the sampled far-field ray.
    assert np.all(cos_align > 0.0)

    # Alignment should improve (or at least not degrade materially) with radius.
    assert np.all(np.diff(cos_align) >= -2.0e-2)

    # Last point should be strongly radial in the far field.
    assert cos_align[-1] > 0.98


def test_bessel_far_field_force_aligns_with_negative_rhat_for_two_disks():
    g = _configure_two_bessel_disks()

    # Mirror the mixed-component scenario and verify attractive far-field
    # direction on an oblique ray.
    radii = np.array([60.0, 120.0, 220.0, 320.0], dtype=float)
    theta = 0.30
    x = radii * np.cos(theta)
    y = np.zeros_like(x)
    z = radii * np.sin(theta)

    ax, ay, az = g.force(x, y, z)
    force = np.column_stack((ax, ay, az))
    pos = np.column_stack((x, y, z))

    minus_rhat = -pos / np.linalg.norm(pos, axis=1, keepdims=True)
    force_norm = np.linalg.norm(force, axis=1)

    assert np.all(force_norm > 0.0)

    cos_align = np.sum(force * minus_rhat, axis=1) / force_norm

    assert np.all(cos_align > 0.0)
    assert np.all(np.diff(cos_align) >= -3.0e-2)
    assert cos_align[-1] > 0.95


def test_bessel_far_field_inferred_mass_curves_are_consistent():
    sigma0 = 1.0
    hR = 4.0
    hZ = 0.8
    g = _configure_single_bessel_disk(sigma0=sigma0, hR=hR, hZ=hZ, nr=192, nz=160, nk=384)

    radii = np.array([30.0, 50.0, 80.0, 120.0, 180.0], dtype=float)
    theta = 0.33
    m_phi, m_a = _far_field_inferred_masses(g, radii, theta)

    m_true = _disk_total_mass(sigma0, hR)

    # Both estimators should remain positive and approach true enclosed mass.
    assert np.all(m_phi > 0.0)
    assert np.all(m_a > 0.0)

    np.testing.assert_allclose(m_phi[-1], m_true, rtol=5.0e-2, atol=1.0e-6)
    np.testing.assert_allclose(m_a[-1], m_true, rtol=5.0e-2, atol=1.0e-6)

    # The two asymptotic mass estimators should agree with one another.
    np.testing.assert_allclose(m_phi[-2:], m_a[-2:], rtol=2.0e-2, atol=1.0e-6)


def test_bessel_interior_poisson_residual_map_is_bounded():
    sigma0 = 1.0
    hR = 4.0
    hZ = 0.8
    g = _configure_single_bessel_disk(sigma0=sigma0, hR=hR, hZ=hZ, nr=192, nz=160, nk=384)

    diag = _poisson_residual_map(g, sigma0, hR, hZ)
    rel = diag["rel"]
    source = np.abs(diag["source"])

    # Ignore an outer belt where truncation/interpolation closure dominates.
    rel_core_grid = rel[3:-3, 3:-3]
    source_core = source[3:-3, 3:-3]

    # Restrict residual statistics to cells where source strength is meaningful.
    source_mask = source_core > 1.0e-3 * np.nanmax(source_core)
    rel_core = rel_core_grid[source_mask]

    assert rel_core.size > 0

    assert np.nanmedian(rel_core) < 9.0e-1
    assert np.nanpercentile(rel_core, 95.0) < 3.5


def test_bessel_matches_direct_thick_disk_reference():
    sigma0 = 1.0
    hR = 4.0
    hZ = 0.8
    g = _configure_single_bessel_disk(sigma0=sigma0, hR=hR, hZ=hZ, nr=192, nz=160, nk=384)

    x = np.array([4.0], dtype=float)
    y = np.array([0.0], dtype=float)
    z = np.array([0.8], dtype=float)

    phi_bessel = g.potential(x, y, z)[0]
    phi_ref = _direct_reference_potential(x[0], z[0], sigma0, hR, hZ)

    np.testing.assert_allclose(phi_bessel, phi_ref, rtol=1.5e-1, atol=2.0e-4)


def test_bessel_resolution_convergence_force_and_potential():
    sigma0 = 1.0
    hR = 4.0
    hZ = 0.8

    x = np.array([2.0, 4.0, 8.0], dtype=float)
    y = np.array([0.0, 0.0, 0.0], dtype=float)
    z = np.array([0.2, 0.8, 1.5], dtype=float)

    g_low = _configure_single_bessel_disk(sigma0=sigma0, hR=hR, hZ=hZ, nr=64, nz=64, nk=128)
    phi_low = g_low.potential(x, y, z)
    ax_low, _, az_low = g_low.force(x, y, z)

    g_mid = _configure_single_bessel_disk(sigma0=sigma0, hR=hR, hZ=hZ, nr=128, nz=128, nk=256)
    phi_mid = g_mid.potential(x, y, z)
    ax_mid, _, az_mid = g_mid.force(x, y, z)

    g_high = _configure_single_bessel_disk(sigma0=sigma0, hR=hR, hZ=hZ, nr=192, nz=160, nk=384)
    phi_high = g_high.potential(x, y, z)
    ax_high, _, az_high = g_high.force(x, y, z)

    err_lm_phi = np.linalg.norm(phi_low - phi_mid)
    err_mh_phi = np.linalg.norm(phi_mid - phi_high)

    force_low = np.column_stack((ax_low, az_low))
    force_mid = np.column_stack((ax_mid, az_mid))
    force_high = np.column_stack((ax_high, az_high))

    err_lm_force = np.linalg.norm(force_low - force_mid)
    err_mh_force = np.linalg.norm(force_mid - force_high)

    # Convergence gate: higher resolution should move less than lower resolution.
    assert err_mh_phi < 0.9 * err_lm_phi
    assert err_mh_force < 0.9 * err_lm_force

    # Absolute sanity: medium and high should already be close for this test set.
    np.testing.assert_allclose(phi_mid, phi_high, rtol=1.5e-1, atol=2.0e-4)
    np.testing.assert_allclose(force_mid, force_high, rtol=2.0e-1, atol=2.0e-4)