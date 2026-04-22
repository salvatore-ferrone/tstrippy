"""
A broad helper module for aiding with construcing orbits and the like
"""

import numpy as np 

def apoapsis_shot(potential, params, radius, eccen, theta, phi, polar_mix, getVecs=False):
    """
    Generate initial position and velocity for an orbit launched from apogee.
    
    Uses pseudo-Keplerian parameters to construct phase-space coordinates
    in the Milky Way potential.
    
    Parameters
    ----------
    potential : function 
        tstrippy.potentials.pouliasis2017pii()
    params : array-like
        Potential parameters (from tstrippy.Parsers.pouliasis2017pii())
    radius : float
        Distance from galactic center (kpc)
    eccen : float
        Eccentricity-like speed scaling, 0 ≤ eccen < 1.
        At eccen=0, speed = circular speed. At eccen→1, speed→0.
    theta : float
        Polar angle in spherical coordinates (radians), 0 ≤ theta ≤ π.
        theta=π/2 is the equatorial plane.
    phi : float
        Azimuthal angle in spherical coordinates (radians), 0 ≤ phi < 2π.
    polar_mix : float
        Polar angle (radians), 0 ≤ polar_mix ≤ π.
        polar_mix=0 gives pure azimuthal (planar) motion.
        polar_mix=π/2 gives equal azimuthal and meridional mix.
        polar_mix=π gives pure meridional motion.
    getVecs : bool, optional
        If True, also return basis vectors. Default False.
    
    Returns
    -------
    rvec : ndarray, shape (3,)
        Position vector [x, y, z]
    vvec : ndarray, shape (3,)
        Velocity vector [vx, vy, vz]
    azimuth_hat : ndarray, shape (3,) (if getVecs=True)
    meridional_hat : ndarray, shape (3,) (if getVecs=True)
    """
    # Build Cartesian position from spherical coordinates
    x = radius * np.cos(phi) * np.sin(theta)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(theta)
    rvec = np.array([x, y, z])
    
    # Compute circular speed from the local force
    fx, fy, fz, _ = potential(params, x, y, z)
    fmag = np.sqrt(fx**2 + fy**2 + fz**2)
    vcirc = np.sqrt(radius * fmag)
    
    # Define spherical basis vectors at (theta, phi)
    azimuth_hat = np.array([-np.sin(phi), np.cos(phi), 0])
    meridional_hat = np.array([np.cos(theta)*np.cos(phi), 
                            np.cos(theta)*np.sin(phi), 
                            -np.sin(theta)])  # Remove the negative sign

    # Mix azimuthal and meridional directions using polar_mix as angle
    vdirection = np.cos(polar_mix) * azimuth_hat + np.sin(polar_mix) * meridional_hat
    vdirection /= np.linalg.norm(vdirection)    
    
    # Scale by eccentricity and circular speed
    vvec = vcirc * (1 - eccen) * vdirection
    
    if getVecs:
        return rvec, vvec, azimuth_hat, meridional_hat
    else:
        return rvec, vvec