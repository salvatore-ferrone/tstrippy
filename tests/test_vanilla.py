"""
Tests the simplest full case 

1. sampling a plummer sphere
2. integrating the host 
3. loading the host's orbit to the integrator
4. integrating the particles positions 

"""

import numpy as np
import tstrippy
import matplotlib as mpl 
import matplotlib.pyplot as plt
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Computer Modern']
mpl.rcParams['text.usetex'] = True
mpl.rcParams['figure.facecolor'] = 'black'
mpl.rcParams['axes.facecolor'] = 'black'
mpl.rcParams['axes.edgecolor'] = 'white'
mpl.rcParams['text.color'] = 'white'
mpl.rcParams['axes.labelcolor'] = 'white'
mpl.rcParams['xtick.color'] = 'white'
mpl.rcParams['ytick.color'] = 'white'
mpl.rcParams['axes.grid'] = True
mpl.rcParams['grid.color'] = 'gray'
mpl.rcParams['grid.alpha'] = 0.3

def test_host_perturber_reinit_smoke():
    G = 1.0

    # Static galaxy: simple Plummer potential
    M_gal = 1e6
    a_gal = 5.0
    staticgalaxy = ["plummer", np.array([G, M_gal, a_gal], dtype=np.float64)]
    # get the crossing time of the galaxy
    tau_gal = np.sqrt ( a_gal ** 3 / (G*M_gal))

    integration_time = 10 * 2 * np.pi * tau_gal

    # Host Plummer properties
    M_host = 1
    a_half_mass = 5e-3
    # get the half mass radius
    a_plummer = tstrippy.code.ergodic.convertHalfMassRadiusToPlummerRadius(a_half_mass)
    # Cluster properties
    N_particles = 100

    # Time setup 
    # the timestep should be based on the cluster's dynamical time 
    tau = np.sqrt(a_plummer**3 / (G*M_host)) 
    dt = 1.0e-2 * tau
    
    nsteps = int(integration_time / dt )
    integrationparameters = [0.0, dt, nsteps]

    # Put the host on an approximately circular orbit in the galaxy Plummer potential
    eccen = 0.9
    R0 = 2.0*a_gal
    x0 = np.array([R0], dtype=np.float64)
    y0 = np.array([0.0], dtype=np.float64)
    z0 = np.array([0.0], dtype=np.float64)

    vc = eccen*np.sqrt(G * M_gal * R0**2 / (R0**2 + a_gal**2)**1.5)
    vx0 = np.array([0.0], dtype=np.float64)
    vy0 = np.array([vc], dtype=np.float64)
    vz0 = np.array([0.0], dtype=np.float64)

    # 1. Find the orbit of the host point in the Plummer potential
    tstrippy.integrator.setstaticgalaxy(*staticgalaxy)
    tstrippy.integrator.setintegrationparameters(*integrationparameters)
    tstrippy.integrator.setinitialkinematics(x0, y0, z0, vx0, vy0, vz0)
    xH, yH, zH, vxH, vyH, vzH = tstrippy.integrator.leapfrogintime(nsteps, 1)
    timeH = tstrippy.integrator.timestamps.copy()

    # 2-3. Save orbit and deallocate
    tstrippy.integrator.deallocate()

    # Collapse host trajectory from shape (1, nsteps+1) to (nsteps+1,)
    xH = xH[0]
    yH = yH[0]
    zH = zH[0]
    vxH = vxH[0]
    vyH = vyH[0]
    vzH = vzH[0]

    # 4. Sample a 10-particle Plummer sphere
    xrel, yrel, zrel, vxrel, vyrel, vzrel = tstrippy.ergodic.isotropicplummer(
        G, M_host, a_half_mass, N_particles
    )

    # 5. Place the cluster at the host initial phase-space point
    xC = xrel + xH[0]
    yC = yrel + yH[0]
    zC = zrel + zH[0]
    vxC = vxrel + vxH[0]
    vyC = vyrel + vyH[0]
    vzC = vzrel + vzH[0]

    # 6. Set up integrator again with host perturber
    tstrippy.integrator.setstaticgalaxy(*staticgalaxy)
    tstrippy.integrator.setintegrationparameters(*integrationparameters)
    tstrippy.integrator.setinitialkinematics(xC, yC, zC, vxC, vyC, vzC)

    tstrippy.integrator.inithostmass("constant", np.array([M_host], dtype=np.float64))
    tstrippy.integrator.inithostradius(a_plummer)
    tstrippy.integrator.inithostkinematics(timeH, xH, yH, zH, vxH, vyH, vzH)

    # 7. Integrate cluster + host together
    tstrippy.integrator.leapfrogtofinalpositions()
    # extract the end points
    x = tstrippy.integrator.xf.copy()
    y = tstrippy.integrator.yf.copy()
    z = tstrippy.integrator.zf.copy()
    vx = tstrippy.integrator.vxf.copy()
    vy = tstrippy.integrator.vyf.copy()
    vz = tstrippy.integrator.vzf.copy()

    # Basic smoke assertions
    assert x.shape[0] == N_particles
    assert y.shape[0] == N_particles
    assert z.shape[0] == N_particles
    assert np.isfinite(x).all()
    assert np.isfinite(y).all()
    assert np.isfinite(z).all()
    assert np.isfinite(vx).all()
    assert np.isfinite(vy).all()
    assert np.isfinite(vz).all()

    tstrippy.integrator.deallocate()


    # make a figure showing the final stream 

    fig, axis = plt.subplots(1, 1, figsize=(8.2, 6.2))
    axis.plot(xH, yH, color="w", linewidth=1)
    axis.scatter(x, y, s=5, color="red")
    # also plot the initial coordinates
    axis.scatter(xC, yC, s=6, color="green")
    factor = 3
    AXIS = {'xlabel': "X [a_gal]", "ylabel": " Y [a_gal]", "aspect": "equal", "xlim": [-factor*a_gal, factor*a_gal], "ylim": [-factor*a_gal, factor*a_gal]}
    axis.set(**AXIS)
    
    # Inset for initial positions
    ax_inset_init = axis.inset_axes([1.05, 0.15, 0.3, 0.3])
    ax_inset_init.scatter(xC, yC, s=6, color="green")
    ax_inset_init.plot(xH[0:1000], yH[0:1000], color="w", linewidth=1)
    inset_lim = 10 * a_plummer
    ax_inset_init.set_xlim(xH[0] - inset_lim, xH[0] + inset_lim)
    ax_inset_init.set_ylim(yH[0] - inset_lim, yH[0] + inset_lim)
    ax_inset_init.set_facecolor('black')
    ax_inset_init.tick_params(colors='white', labelsize=8)
    
    # Inset for final positions
    ax_inset_final = axis.inset_axes([1.01, 0.55, 0.3, 0.3])
    ax_inset_final.scatter(x, y, s=5, color="red")
    ax_inset_final.plot(xH[-1000:-1], yH[-1000:-1],color="w",linewidth=1)
    ax_inset_final.set_xlim(xH[-1] - inset_lim, xH[-1] + inset_lim)
    ax_inset_final.set_ylim(yH[-1] - inset_lim, yH[-1] + inset_lim)
    ax_inset_final.set_facecolor('black')
    ax_inset_final.tick_params(colors='white', labelsize=8)
    
    fig.savefig("../test2.png")



    return None