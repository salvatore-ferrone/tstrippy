import tstrippy
from astropy import units as u
from astropy import constants as const
from astropy import coordinates as coord
import numpy as np

absolutetolerance=1e-10



def test_vanilla_clusters():
    MWparams        = tstrippy.Parsers.potential_parameters.pouliasis2017pii()
    MWrefframe      = tstrippy.Parsers.potential_parameters.MWreferenceframe()  
    x,y,z,vx,vy,vz  = load_globular_clusters_in_galactic_coordinates(MWrefframe)
    staticgalaxy    = ["pouliasis2017pii", MWparams]
    integrationtime =   5e9 * u.yr
    timestep        =   1e5 * u.yr
    initialkinematics=[x,y,z,vx,vy,vz]
    backwardOrbit,forwardOrbit=\
        vanilla_clusters(integrationtime,timestep,staticgalaxy,initialkinematics)
    assert len(backwardOrbit)==7
    assert len(forwardOrbit)==7
    dr = np.sqrt((backwardOrbit[1]-forwardOrbit[1])**2 +\
                    (backwardOrbit[2]-forwardOrbit[2])**2 +\
                    (backwardOrbit[3]-forwardOrbit[3])**2)
    dv = np.sqrt((backwardOrbit[4]-forwardOrbit[4])**2 +\
                    (backwardOrbit[5]-forwardOrbit[5])**2 +\
                    (backwardOrbit[6]-forwardOrbit[6])**2)
    
    assert np.all(dr<absolutetolerance)
    assert np.all(dv<absolutetolerance)
    return None


def vanilla_clusters(integrationtime,timestep,staticgalaxy,initialkinematics):
    """
    do the backward and forward integration of the vanilla clusters
    """
    assert isinstance(integrationtime,u.Quantity)
    assert isinstance(timestep,u.Quantity)
    unitT, unitV, unitD, unitM, unitG, G = loadunits()
    Ntimestep=int(integrationtime.value/timestep.value)
    dt=timestep.to(unitT)
    currenttime=0*unitT
    integrationparameters=[currenttime.value,dt.value,Ntimestep]
    nObj = initialkinematics[0].shape[0]

    tstrippy.integrator.setstaticgalaxy(*staticgalaxy)
    tstrippy.integrator.setinitialkinematics(*initialkinematics)
    tstrippy.integrator.setintegrationparameters(*integrationparameters)
    tstrippy.integrator.setbackwardorbit()
    xBackward,yBackward,zBackward,vxBackward,vyBackward,vzBackward=\
        tstrippy.integrator.leapfrogintime(Ntimestep,nObj)
    tBackward=tstrippy.integrator.timestamps.copy()
    tstrippy.integrator.deallocate()

    #### Now compute the orbit forward
    currenttime=tBackward[-1]*unitT

    integrationparameters=[currenttime.value,dt.value,Ntimestep]
    x0,y0,z0=xBackward[:,-1],yBackward[:,-1],zBackward[:,-1]
    vx0,vy0,vz0 = -vxBackward[:,-1],-vyBackward[:,-1],-vzBackward[:,-1]
    initialkinematics=[x0,y0,z0,vx0,vy0,vz0]
    tstrippy.integrator.setstaticgalaxy(*staticgalaxy)
    tstrippy.integrator.setintegrationparameters(*integrationparameters)
    tstrippy.integrator.setinitialkinematics(*initialkinematics)
    xForward,yForward,zForward,vxForward,vyForward,vzForward=\
        tstrippy.integrator.leapfrogintime(Ntimestep,nObj)
    tForward=tstrippy.integrator.timestamps.copy()
    tstrippy.integrator.deallocate()

    # flip the backorbits such that the point in the past, 
    # which should be the common starting point, 
    # is the first point for both the forward and backward orbits
    tBackward=tBackward[::-1]
    xBackward,yBackward,zBackward=xBackward[:,::-1],yBackward[:,::-1],zBackward[:,::-1]
    vxBackward,vyBackward,vzBackward=-vxBackward[:,::-1],-vyBackward[:,::-1],-vzBackward[:,::-1]
    backwardOrbit  = [tBackward,xBackward,yBackward,zBackward,vxBackward,vyBackward,vzBackward]
    forwardOrbit   = [tForward,xForward,yForward,zForward,vxForward,vyForward,vzForward]
    return backwardOrbit,forwardOrbit



    

def load_globular_clusters_in_galactic_coordinates(MWrefframe):
    unitT, unitV, unitD, unitM, unitG, G = loadunits()
    GCdata  =   tstrippy.Parsers.baumgardtMWGCs().data
    skycoordinates=coord.SkyCoord(
        ra=GCdata['RA'],
        dec=GCdata['DEC'],
        distance=GCdata['Rsun'],
        pm_ra_cosdec=GCdata['mualpha'],
        pm_dec=GCdata['mu_delta'],
        radial_velocity=GCdata['RV'],)
    galacticcoordinates = skycoordinates.transform_to(MWrefframe)
    x,y,z=galacticcoordinates.cartesian.xyz.to(unitD).value
    vx,vy,vz=galacticcoordinates.velocity.d_xyz.to(unitV).value
    return x,y,z,vx,vy,vz


def bar_movement_ferrone2023():
    """
    Given the reference frame of the Galaxy, the bar OMEGA is negative 
    """

    # oreitnation and bar pattern speed
    unitT, unitV, unitD, unitM, unitG, G=loadunits()
    theta0= 25 * (np.pi/180) 
    omega =  28  * 2*np.pi * unitV / unitD
    omega = -omega.value
    barpolycoeff=[theta0,omega]
    return barpolycoeff

def barparams_ferrone2023():
    unitT, unitV, unitD, unitM, unitG, G=loadunits()
    Mbar = 990.0*2.32*1e7 * unitM    
    abar = 4 * unitD
    bbar = 1 * unitD
    cbar = 0.5 * unitD
    Mbar = 990.0*2.32*1e7 * unitM
    barparams = [G,Mbar.value,abar.value,bbar.value,cbar.value]    
    return barparams



def loadunits():
    # Load the units
    unitbasis = tstrippy.Parsers.potential_parameters.unitbasis
    unitT=u.Unit(unitbasis['time'])
    unitV=u.Unit(unitbasis['velocity'])
    unitD=u.Unit(unitbasis['distance'])
    unitM=u.Unit(unitbasis['mass'])
    unitG=u.Unit(unitbasis['G'])
    G = const.G.to(unitG).value
    return unitT, unitV, unitD, unitM, unitG, G
