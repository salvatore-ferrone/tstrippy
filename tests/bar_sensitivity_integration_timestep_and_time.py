import tstrippy
import astropy.units as u
import astropy.constants as const
import astropy.coordinates as coord
import numpy as np
import matplotlib.pyplot as plt
import datetime



def main(targetGC):
    # Load the units
    unitT, unitV, unitD, unitM, unitG, G=loadunits()
    # Load the galaxy parameters
    MWparams, MWrefframe = loadGalaxy()
    # Load the bar parameters
    barparams = barparams_ferrone2023()
    # Load the bar movement
    barpolycoeff = bar_movement_ferrone2023()
    # Reduce the mass of the disks to account for the bar
    MWparams[5] = 1120.0 * 2.32*10**7 
    MWparams[8] = 1190.0 * 2.32*10**7 

    x0,y0,z0,vx0,vy0,vz0 = pick_globular_cluster(targetGC, MWrefframe)
    Galaxy = ["pouliasis2017pii", MWparams]
    bar = ["longmuralibar", barparams, barpolycoeff]
    initialkinematics = [x0,y0,z0,vx0,vy0,vz0]
    # Pick the target globular cluster
    
    integrationtime = 1e9
    timestep = 1e7
    backward_and_forward_orbit(integrationtime,timestep)
    
    return None


def backward_and_forward_orbit(integartionparams,integrationparameters,staticgalaxy,initialkinematics,galacticbar):
    
    integrationtime,dt,Ntimestep=integrationparameters

    assert integrationTime > 0
    assert timestep > 0
    assert isinstance(integrationTime, (int, float))
    assert isinstance(timestep, (int, float))
    unitT, unitV, unitD, unitM, unitG, G = loadunits()

    Ntimestep=int(integrationTime/timestep)
    T,dt=integrationTime*u.yr,timestep*u.yr
    T,dt=T.to(unitT),dt.to(unitT)
    currenttime = 0*unitT

    xBackward,yBackward,zBackward,vxBackward,vyBackward,vzBackward = \
        backward_orbit([currenttime,dt,Ntimestep],MWparams,[x0,y0,z0,vx0,vy0,vz0],[barname,barparams,barpolycoeff])


def backward_orbit(integrationparameters,staticgalaxy,initialkinematics,galacticbar):
    nObj = 1 # only integrating one object
    currenttime,dt,Ntimestep=integrationparameters
    MWname, MWparams = staticgalaxy
    x0,y0,z0,vx0,vy0,vz0 = initialkinematics
    barname,barparams,barpolycoeff = galacticbar

    tstrippy.integrator.setstaticgalaxy(MWname,MWparams)
    tstrippy.integrator.setbackwardorbit()
    tstrippy.integrator.setintegrationparameters(currenttime.value,dt.value,Ntimestep)
    tstrippy.integrator.setinitialconditions(x0,y0,z0,vx0,vy0,vz0)
    tstrippy.integrator.setgalacticbar(barname,barparams,barpolycoeff)
    starttime = datetime.datetime.now()
    xBackward,yBackward,zBackward,vxBackward,vyBackward,vzBackward=tstrippy.integrator.leapfrogintime(Ntimestep,nObj)
    endtime=datetime.datetime.now()
    print("Backward orbit took",endtime-starttime)
    tstrippy.integrator.deallocate()
    return xBackward,yBackward,zBackward,vxBackward,vyBackward,vzBackward



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



def loadGalaxy():
    MWparams = tstrippy.Parsers.potential_parameters.pouliasis2017pii()
    MWrefframe= tstrippy.Parsers.potential_parameters.MWreferenceframe()    
    return MWparams, MWrefframe


def barparams_ferrone2023():
    unitT, unitV, unitD, unitM, unitG, G=loadunits()
    Mbar = 990.0*2.32*1e7 * unitM    
    abar = 4 * unitD
    bbar = 1 * unitD
    cbar = 0.5 * unitD
    Mbar = 990.0*2.32*1e7 * unitM
    barparams = [G,Mbar.value,abar.value,bbar.value,cbar.value]    
    return barparams


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



def pick_globular_cluster(targetGC, MWrefframe):
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

    # Pick the target globular cluster
    clusterindex = np.where(GCdata['Cluster']==targetGC)[0][0]
    x0,y0,z0=x[clusterindex],y[clusterindex],z[clusterindex]
    vx0,vy0,vz0=vx[clusterindex],vy[clusterindex],vz[clusterindex]
    return x0,y0,z0,vx0,vy0,vz0





if __name__ == "__main__":
    targetGC = "NGC5139"
    main(targetGC)