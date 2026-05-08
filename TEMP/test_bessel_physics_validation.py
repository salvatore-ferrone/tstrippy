import math

import numpy as np
import pytest

import gravitymini


def _configure_single_bessel_disk(sigma0=1.0, hR=4.0, hZ=0.8, nr=128, nz=128, nk=256):
    g = gravitymini.gravity
    b = gravitymini.besselbfe

    g.cleargravity()
    b.bessel_init(nr, nz, nk, hR, hZ)
    g.addgravitycomponent("exponentialdisk", [sigma0, hR, hZ])
    g.finalizegravity()
    return g


def _disk_total_mass(sigma0, hR):
    return 2.0 * math.pi * sigma0 * hR**2


def _direct_reference_potential(R_eval, z_eval, sigma0, hR, hZ, nR=140, nZ=120, nphi=180):
    gconst = float(gravitymini.gravity.gravity_g_default)

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
    gconst = float(gravitymini.gravity.gravity_g_default)
    radius = x[0]
    phi_ref = -gconst * mass / radius
    ax_ref = -gconst * mass / radius**2

    np.testing.assert_allclose(phi, phi_ref, rtol=2.5e-1, atol=1e-6)
    np.testing.assert_allclose(ax[0], ax_ref, rtol=2.5e-1, atol=1e-6)
    np.testing.assert_allclose(ay[0], 0.0, atol=1e-12)
    np.testing.assert_allclose(az[0], 0.0, atol=1e-12)


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