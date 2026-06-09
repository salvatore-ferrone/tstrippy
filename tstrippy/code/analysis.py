"""
A module containing functions used for routine analysis
"""

import numpy as np 

def compute_hostcluster_mass_evolution(simulator, present_day_host_mass):
    r"""

    Define the mass evolution as:

    $$ 
    M_p = \frac{M_{\rm{today}}}{N_{\rm{bound}}+1}\\
    $$
    
    Adding one in the denominator is for the center of mass

    $$
    M(t) = M_p \left[ N_{\rm{particles}} + 1 - \int_{t_{0}}^{t} \sum_{i}^{N_{\rm{esc}}}\delta\left(t-t_{\rm{esc},i}\right) dt\right]
    $$

    Note that:

    $$ N_{\rm{particles}} = N_{\rm{bound}} + N_{\rm{esc}}$$

    Limit testing. If all escape:

    $$ M(t_0) = ( N_{\rm{particles}} +1 ) * M_{\rm{today}} $$

    If none escape, then

    $$ M(t_0) =  M_{\rm{today}} $$

    Note that the lists are appended with estimate of the "initial mass" and with the present day mass
    This avoids computation errors in the convergence analysis in the evend that no particles are ionized
    
    """
    # type/interface check
    if simulator.__class__.__name__ != "fortran":
        raise TypeError(f"Expected fortran simulator object, got {simulator.__class__}")

    # value check
    m = float(present_day_host_mass)
    if not np.isfinite(m) or m <= 0.0:
        raise ValueError(f"present_day_host_mass must be finite and > 0, got {present_day_host_mass}")
    

    ionized, escapetime=simulator.get_hostcluster_ionization_state(simulator.nparticles)
    ionized=ionized.astype(np.bool_)
    Nparticles = ionized.shape[0]
    Nbound = (~ionized).sum()
    Nionized = ionized.sum()

    mass_per_particle = present_day_host_mass/(Nbound+1) # add one for the center of mass of the cluster

    host_mass_at_start_of_simulation = mass_per_particle * (Nparticles+1) 

    nescape = np.cumsum(np.ones(Nionized))

    mass_evolution = host_mass_at_start_of_simulation - mass_per_particle * nescape

    mass_evolution_timestamps = np.sort(escapetime[ionized])

    # append both arrays for the initial and final mass 
    finaltime=simulator.timestamps[-1]
    initialtime = simulator.timestamps[0]

    mass_evolution_timestamps = np.concatenate(([initialtime], mass_evolution_timestamps, [finaltime]))
    mass_evolution = np.concatenate(([host_mass_at_start_of_simulation], mass_evolution, [present_day_host_mass]))
    
    return mass_evolution_timestamps, mass_evolution


def transform_galactocentric_to_tail_coordinates(hostorbit_timestamps,hostorbit,stream, t0=0, returnindexes=False):
    """Transform stream phase-space coordinates into local tail coordinates.

    Parameters
    ----------
    hostorbit_timestamps : array-like, shape (ntimestamps,)
        Time values associated with the host orbit samples.
    hostorbit : array-like, shape (6, ntimestamps)
        Host orbit phase-space sampled at ``hostorbit_timestamps``, ordered as
        ``(x, y, z, vx, vy, vz)``.
    stream : array-like, shape (6, nparticles)
        Stream particle phase-space, ordered as ``(x, y, z, vx, vy, vz)``.
    t0 : float, optional
        Reference time used to define zero path length along the orbit.
    returnindexes : bool, optional
        If True, also return nearest-orbit indices for each particle.

    Returns
    -------
    tailcoordinates : tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        Tuple ``(xprime, yprime, zprime, vxprime, vyprime, vzprime)`` with one
        value per particle.
    indexes : np.ndarray, optional
        Returned only when ``returnindexes=True``. Nearest host-orbit sample
        index for each particle.

    Raises
    ------
    ValueError
        If array dimensionality or required shapes are invalid.
    """
    hostorbit_timestamps = np.asarray(hostorbit_timestamps)
    hostorbit = np.asarray(hostorbit)
    stream = np.asarray(stream)

    if hostorbit_timestamps.ndim != 1:
        raise ValueError("hostorbit_timestamps must be a 1D array")

    if hostorbit.ndim != 2:
        raise ValueError("hostorbit must be a 2D array with shape (6, ntimestamps)")
    if hostorbit.shape[0] != 6:
        raise ValueError(f"hostorbit must have shape (6, ntimestamps); got {hostorbit.shape}")
    if hostorbit.shape[1] != hostorbit_timestamps.shape[0]:
        raise ValueError(
            "hostorbit_timestamps length must match hostorbit columns "
            f"({hostorbit_timestamps.shape[0]} != {hostorbit.shape[1]})"
        )

    if stream.ndim != 2:
        raise ValueError("stream must be a 2D array")
    if stream.shape[0] != 6:
        raise ValueError(f"stream must have shape (6, nparticles); got {stream.shape}")

    xORB,yORB,zORB,vxORB,vyORB,vzORB = hostorbit 
    xp,yp,zp,vxp,vyp,vzp = stream 

    indexes=_find_closest_orbital_indices(xORB,yORB,zORB,xp,yp,zp)
    xprimeOrbi=_compute_orbital_path_length(hostorbit_timestamps,xORB,yORB,zORB,t0=t0)

    xprimeP=np.zeros(xp.shape)
    yprimeP=np.zeros(xp.shape)
    zprimeP=np.zeros(xp.shape)
    vxprimeP=np.zeros(xp.shape)
    vyprimeP=np.zeros(xp.shape)
    vzprimeP=np.zeros(xp.shape)
    for i in range(xp.shape[0]):
        POSpart=np.array([xp[i],yp[i],zp[i]])
        POS_near_point=np.array([xORB[indexes[i]],yORB[indexes[i]],zORB[indexes[i]]])
        vel_near= np.array([vxORB[indexes[i]],vyORB[indexes[i]],vzORB[indexes[i]]])
        VelPart = np.array([vxp[i], vyp[i], vzp[i]])
        dV=VelPart-vel_near
        dX = POSpart-POS_near_point
        vunit,runit,zunit,compensation=_get_local_unit_vectors_and_particle_offset_compensation(
            POSpart,POS_near_point,vel_near)
        # outputs
        xprimeP[i]=xprimeOrbi[indexes[i]]+np.dot(compensation,vunit)
        yprimeP[i]=np.dot(dX,runit)
        zprimeP[i]=np.dot(dX,zunit)
        vxprimeP[i] = np.dot(dV,vunit)
        vyprimeP[i] = np.dot(dV,runit)
        vzprimeP[i] = np.dot(dV,zunit)    

    tailcoordinates = (xprimeP, yprimeP, zprimeP, vxprimeP, vyprimeP, vzprimeP)

    if returnindexes:
        return tailcoordinates, indexes

    return tailcoordinates

def window_orbit_by_local_dynamical_time(
    hostorbit_timestamps: np.ndarray,
    hostorbit: np.ndarray,
    time_of_interest: float,
    nDynTimes: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Window orbit coordinates around a time of interest using local dynamical time.
    This is useful since I am doing a naive computation of finding the nearest point of each particle to the orbit.

    Args:
        hostorbit_timestamps (np.ndarray): Array of time values with shape
            (ntimestamps,).
        hostorbit (np.ndarray): Phase-space orbit array with shape
            (6, ntimestamps), ordered as (x, y, z, vx, vy, vz).
        time_of_interest (float): Current time value.
        nDynTimes (int): Number of dynamical times.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Filtered timestamps and filtered host
        orbit array with shape (6, nfiltered).
    """
    
    hostorbit_timestamps = np.asarray(hostorbit_timestamps)
    hostorbit = np.asarray(hostorbit)

    if hostorbit_timestamps.ndim != 1:
        raise ValueError("hostorbit_timestamps must be a 1D array")
    if hostorbit.ndim != 2:
        raise ValueError("hostorbit must be a 2D array with shape (6, ntimestamps)")
    if hostorbit.shape[0] != 6:
        raise ValueError(f"hostorbit must have shape (6, ntimestamps); got {hostorbit.shape}")
    if hostorbit.shape[1] != hostorbit_timestamps.shape[0]:
        raise ValueError(
            "hostorbit_timestamps length must match hostorbit columns "
            f"({hostorbit_timestamps.shape[0]} != {hostorbit.shape[1]})"
        )

    tORB = hostorbit_timestamps
    xtORB, ytORB, ztORB, vxtORB, vytORB, vztORB = hostorbit

    assert ( nDynTimes > 0), "nDynTimes must be greater than 0"
    assert ( nDynTimes < 5), "nDynTimes must be less than 5"
    assert( time_of_interest > np.min(tORB)), "time_of_interest must be greater than the minimum time in tORB"
    assert( time_of_interest < np.max(tORB)), "time_of_interest must be less than the maximum time in tORB"
    
    rORB = np.sqrt(xtORB ** 2 + ytORB ** 2 + ztORB ** 2)
    vORB = np.sqrt(vxtORB ** 2 + vytORB ** 2 + vztORB ** 2)

    todayindex = np.argmin(np.abs(tORB - time_of_interest))
    ttoday = tORB[todayindex]

    Tdyn = np.median(rORB / vORB)

    cond1 = (tORB - ttoday) < (nDynTimes * Tdyn)
    cond2 = (tORB - ttoday) > (-nDynTimes * Tdyn)
    cond = cond1 * cond2

    XORB = xtORB[cond]
    YORB = ytORB[cond]
    ZORB = ztORB[cond]
    VXORB = vxtORB[cond]
    VYORB = vytORB[cond]
    VZORB = vztORB[cond]
    TORB = tORB[cond]

    filtered_hostorbit = np.vstack([XORB, YORB, ZORB, VXORB, VYORB, VZORB])
    return TORB, filtered_hostorbit


def _compute_orbital_path_length(tORB:np.ndarray,xORB:np.ndarray,yORB:np.ndarray,zORB:np.ndarray,t0:float=0):
    """
    Calculate the xprime orbit coordinates.
    
    xprime is the pathlength along the orbit ahead of behind the center of mass.
    ORB refers to the orbit of the globular cluster.

    Parameters:
    - tORB (np.ndarray): Array of time values.
    - xORB (np.ndarray): Array of x-coordinate values.
    - yORB (np.ndarray): Array of y-coordinate values.
    - zORB (np.ndarray): Array of z-coordinate values.
    - t0 (float): Reference time value (default: 0).

    Returns:
    - orbital_path_length (np.ndarray): Array of xprime orbit coordinates.
    """
    
    orbital_path_length=np.zeros(xORB.shape)
    center_of_mass_index=np.argmin(np.abs(tORB-t0))
    dx = np.diff(xORB-xORB[center_of_mass_index])
    dy = np.diff(yORB-yORB[center_of_mass_index])
    dz = np.diff(zORB-zORB[center_of_mass_index])
    dr = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    ahead = dr[center_of_mass_index::]
    behind = dr[:center_of_mass_index]
    orbital_path_length[center_of_mass_index+1::] = np.cumsum(ahead)
    orbital_path_length[:center_of_mass_index] = -np.flip(np.cumsum(np.flip(behind)))
    return orbital_path_length

def _find_closest_orbital_indices(xORB:np.ndarray,yORB:np.ndarray,zORB:np.ndarray,xp:np.ndarray,yp:np.ndarray,zp:np.ndarray) -> np.ndarray:
    """
    Returns an array of indexes representing the closest orbital index for each point in the given arrays.

    Parameters:
    - xORB (np.ndarray): Array of x-coordinates of orbital points.
    - yORB (np.ndarray): Array of y-coordinates of orbital points.
    - zORB (np.ndarray): Array of z-coordinates of orbital points.
    - xp (np.ndarray): Array of x-coordinates of points to find closest orbital index for.
    - yp (np.ndarray): Array of y-coordinates of points to find closest orbital index for.
    - zp (np.ndarray): Array of z-coordinates of points to find closest orbital index for.

    Returns:
    - indexes (np.ndarray): Array of nearest corresponding orbital index for each particle.


    For instance, the nearst time stamp along the orbit to a given particle is:
    xp[i], yp[i], zp[i] is the particle position
    tORB[indexes[i]], xORB[indexes[i]], yORB[indexes[i]], zORB[indexes[i]]
    the length of the indexes is the same as the particles

    """
    indexes=np.zeros(xp.shape,dtype=int)
    for i in range(xp.shape[0]):
        dx,dy,dz = xp[i]-xORB,yp[i]-yORB,zp[i]-zORB
        dist=dx**2 + dy**2 + dz**2
        indexes[i]=np.argmin(dist) 
    return indexes

def _get_local_unit_vectors_and_particle_offset_compensation(particle_position: np.ndarray,position_orbital_near_point: np.ndarray,velocity_orbital_near_point: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    '''
    def get_local_unit_vectors_and_particle_offset_compensation(particle_position, position_orbital_near_point, velocity_orbital_near_point):
        """
        The coordinate system of tail coordinates is defined locally for each particle. 
        This function finds those three unit vectors. They are defined as 
        runit: unit radial vector from the galactic center to the nearest point along the orbit to the particle
        vunit: unit velocity vector in the direction of the orbital motion of the host globular cluster at the nearest point along the orbit to the particle
        zunit: unit vector perpendicular to the plane defined by runit and vunit. The angular momentum in essence.

        Parameters:
        - particle_position (np.ndarray): Vector with 3 components representing the particle position.
        - position_orbital_near_point (np.ndarray): Vector with 3 components representing the position of the orbital near point.
        - velocity_orbital_near_point (np.ndarray): Vector with 3 components representing the velocity of the orbital near point.

        Returns:
        - vunit (np.ndarray): Vector representing the unit velocity vector.
        - runit (np.ndarray): Vector representing the unit radial vector.
        - zunit (np.ndarray): Vector representing the unit z vector.
        - compensation (np.ndarray): Vector representing the particle offset compensation.
        """
        # Implementation goes here
    '''
    assert particle_position.shape[0]==3,"particle_position must be vector with 3 components"
    assert position_orbital_near_point.shape[0]==3,"position_orbital_near_point must be vector with 3 components"
    assert velocity_orbital_near_point.shape[0]==3,"velocity_orbital_near_point must be vector with 3 components"
    dX = particle_position-position_orbital_near_point
    vunit =velocity_orbital_near_point/np.linalg.norm(velocity_orbital_near_point)
    # adjust the point slightly, due to the discrete grid
    compensation=np.dot(dX,vunit)*vunit
    interpolationPoint=position_orbital_near_point+compensation 
    runit = interpolationPoint/np.linalg.norm(interpolationPoint)
    zunit=np.cross(runit,vunit)
    return vunit,runit,zunit,compensation
 
