"""
tstrippy/code/bfe.py

Basis Function Expansion (BFE) utilities for axisymmetric potentials.

These functions support working with the Legendre / spherical-harmonic
basis expansion implemented in the Fortran ``potentials`` module.  They
provide NumPy-only post-processing (density reconstruction, Legendre
evaluation) as well as thin wrappers around the Fortran initialisation
and evaluation routines that handle array layout and dtype conversion.
"""

import numpy as np


def legendre_all(mu, lmax):
    """Evaluate Legendre polynomials P_l(mu) for l = 0 … lmax.

    Uses the three-term recurrence relation, matching the convention
    inside the Fortran ``mathutils`` module.

    Parameters
    ----------
    mu : array_like, shape (N,)
        Cosine of the polar angle, flattened to 1-D.
    lmax : int
        Maximum degree.

    Returns
    -------
    p : ndarray, shape (lmax+1, N)
        ``p[l, i]`` is P_l(mu[i]).
    """
    mu = np.asarray(mu, dtype=float).ravel()
    p = np.zeros((lmax + 1, mu.size), dtype=float)
    p[0] = 1.0
    if lmax == 0:
        return p
    p[1] = mu
    for ell in range(1, lmax):
        p[ell + 1] = ((2 * ell + 1) * mu * p[ell] - ell * p[ell - 1]) / (ell + 1)
    return p


def reconstruct_rho_from_basis(x, y, z, basis_r_grid, basis_rho_l_grid, lmax):
    """Reconstruct the density field from the Legendre basis tables.

    Interpolates each rho_l(r) coefficient in log-r space and sums
    the even-l Legendre series::

        rho(r, mu) = sum_{l=0,2,4,...}^{lmax} rho_l(r) * P_l(mu)

    Parameters
    ----------
    x, y, z : array_like
        Cartesian coordinates (arbitrary broadcastable shape).
    basis_r_grid : ndarray, shape (nr,)
        Radial grid used in the basis expansion (must be log-spaced).
    basis_rho_l_grid : ndarray, shape (lmax+1, nr)
        Legendre coefficients rho_l(r) at each grid point, as stored
        in ``tstrippy.potentials.basis_rho_l_grid``.
    lmax : int
        Maximum Legendre degree used in the expansion.

    Returns
    -------
    rho : ndarray
        Reconstructed density, same shape as the input arrays.
    """
    shape = np.shape(x)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)

    r = np.sqrt(x * x + y * y + z * z)
    r_clamped = np.clip(r, basis_r_grid[0], basis_r_grid[-1])
    mu = np.zeros_like(r_clamped)
    mask = r_clamped > 0.0
    mu[mask] = z[mask] / r_clamped[mask]
    mu = np.clip(mu, -1.0, 1.0)

    rf = r_clamped.ravel()
    muf = mu.ravel()
    lr = np.log(basis_r_grid)
    lrf = np.log(rf)

    p = legendre_all(muf, lmax)
    rho = np.zeros_like(rf)
    for ell in range(0, lmax + 1, 2):
        rho_l = np.interp(lrf, lr, basis_rho_l_grid[ell, :])
        rho += rho_l * p[ell]
    return rho.reshape(shape)


def init_basis_tables(potential_fn, params, G, lmax, nr, rmin, rmax):
    """Initialise the Fortran basis tables for an axisymmetric potential.

    Clears any existing expansion, initialises a log-spaced radial grid,
    then triggers the density projection and potential table computation
    via a single dummy evaluation of the supplied potential function.

    Parameters
    ----------
    potential_fn : callable
        The Fortran potential subroutine. It must accept
        ``(params, x, y, z)`` and return ``(ax, ay, az, phi)``.
    params : array_like, shape (N,)
        Parameter array passed directly to ``potential_fn``.
    G : float
        Gravitational constant in the user's unit system.
    lmax : int
        Maximum Legendre degree for the expansion.
    nr : int
        Number of radial grid points.
    rmin, rmax : float
        Minimum and maximum radii for the log-spaced grid [kpc].

    Returns
    -------
    basis_r_grid : ndarray, shape (nr,)
        Radial grid stored in the Fortran module after initialisation.
    basis_rho_l_grid : ndarray, shape (lmax+1, nr)
        Legendre density coefficients after projection.
    """
    import tstrippy

    r_grid = np.logspace(np.log10(rmin), np.log10(rmax), nr)
    tstrippy.potentials.clearaxisymmetricbasisexpansion()
    tstrippy.potentials.initaxisymmetricbasisexpansion(G, int(lmax), r_grid)

    x0 = np.array([1.0], dtype=float)
    y0 = np.array([0.0], dtype=float)
    z0 = np.array([0.0], dtype=float)
    potential_fn(params, x0, y0, z0)

    basis_r_grid = np.asarray(tstrippy.potentials.basis_r_grid, dtype=float)
    basis_rho_l_grid = np.asarray(tstrippy.potentials.basis_rho_l_grid, dtype=float)
    return basis_r_grid, basis_rho_l_grid


def evaluate_phi_map(potential_fn, params, x, y, z):
    """Evaluate the gravitational potential on an arbitrary coordinate array.

    Handles flattening, dtype coercion, and reshaping so callers can pass
    2-D meshgrids directly to any Fortran potential subroutine.

    Parameters
    ----------
    potential_fn : callable
        The Fortran potential subroutine, e.g.
        ``tstrippy.potentials.exponential_oblate_halo``.  Must accept
        ``(params, x, y, z)`` and return ``(ax, ay, az, phi)``.
    params : array_like
        Parameter array passed directly to ``potential_fn``.
    x, y, z : array_like
        Cartesian coordinates (arbitrary broadcastable shape).

    Returns
    -------
    phi : ndarray
        Gravitational potential, same shape as the input arrays.
    """
    shape = np.shape(x)
    xf = np.asarray(x, dtype=float).ravel()
    yf = np.asarray(y, dtype=float).ravel()
    zf = np.asarray(z, dtype=float).ravel()
    _, _, _, phi = potential_fn(params, xf, yf, zf)
    return np.asarray(phi, dtype=float).reshape(shape)
