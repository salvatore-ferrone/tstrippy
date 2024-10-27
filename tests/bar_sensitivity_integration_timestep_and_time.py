import tstrippy
import astropy.units as u
import astropy.constants as const
import astropy.coordinates as coord
import numpy as np
import matplotlib.pyplot as plt
import datetime
import matplotlib.colors as mcolors



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
    # Extract the initial conditions of the target globular cluster
    x0,y0,z0,vx0,vy0,vz0 = pick_globular_cluster(targetGC, MWrefframe)

    #### Begin storing the input data for the integrator
    Galaxy = ["pouliasis2017pii", MWparams]
    bar = ["longmuralibar", barparams, barpolycoeff]
    initialkinematics = [x0,y0,z0,vx0,vy0,vz0]
    # Pick the target globular cluster
    integrationtime = 1e9
    timestep = 1e7
    integrationparameters = [integrationtime,timestep]
    backward_and_forward_orbit(integrationtime,timestep,Galaxy,initialkinematics,bar)
    
    return None


def backward_and_forward_orbit(integrationtime,timestep,staticgalaxy,initialkinematics,galacticbar):
    

    assert integrationtime > 0
    assert timestep > 0
    assert isinstance(integrationtime, (int, float))
    assert isinstance(timestep, (int, float))
    unitT, unitV, unitD, unitM, unitG, G = loadunits()

    # convert to integration units
    Ntimestep=int(integrationtime/timestep)
    T,dt=integrationtime*u.yr,timestep*u.yr
    T,dt=T.to(unitT),dt.to(unitT)
    # make sure the current time is set to today for backward integration
    currenttime = 0*unitT
    integrationparameters = [currenttime,dt,Ntimestep]
    xBackward,yBackward,zBackward,vxBackward,vyBackward,vzBackward = \
        backward_orbit(integrationparameters,staticgalaxy,initialkinematics,galacticbar)
    

    print("len(xBackward)",len(xBackward),)
    plotdata0={"x":xBackward[0],"y":yBackward[0]}
    pltobj0={"label":"Backward orbit","color":"blue"}

    plot_orbits([plotdata0],[pltobj0])
    return None


def plot_orbits(plotData,plotProperties,axisconfig={
        'xlabel':"X [kpc]",
        'ylabel':"Y [kpc]",
        'aspect':'equal',},
        legendconfig={'loc':'center left', 'bbox_to_anchor':(1.05, 0.5)},
        ):

    assert isinstance(plotData, list)
    assert isinstance(plotProperties, list)
    assert isinstance(axisconfig, dict)
    assert len(plotData) == len(plotProperties)

    nplots = len(plotData)
    scattersizes = np.logspace(2,0,nplots)
    linewidths = np.linspace(4,0.5,nplots)
    fig,axis=plt.subplots(1,1,figsize=(6,5))
    for i in range(nplots):
        plotProperties[i]['linewidth']=linewidths[i]
        plotProperties[i]['zorder']=nplots-i
        axis.plot(plotData[i]['x'],plotData[i]['y'],**plotProperties[i])


    cmap = plt.get_cmap('tab10') 
    for i in range(nplots):
        mylabel = plotProperties[i].get('label','')
        axis.scatter(plotData[i]['x'][0],plotData[i]['y'][0],s=scattersizes[i],label=mylabel+" Initial",color=plotProperties[i]['color'])
        mycolor = cmap(i % cmap.N)
        axis.scatter(plotData[i]['x'][-1],plotData[i]['y'][-1],s=scattersizes[i],label=mylabel+" Final",color=mycolor)

    axis.legend(**legendconfig)
    axis.set(**axisconfig)
    fig.tight_layout()
    fig.savefig("orbit_integration_timestep_and_time.png")



def backward_orbit(integrationparameters,staticgalaxy,initialkinematics,galacticbar):
    nObj = 1 # only integrating one object
    currenttime,dt,Ntimestep=integrationparameters
    MWname, MWparams = staticgalaxy
    x0,y0,z0,vx0,vy0,vz0 = initialkinematics
    barname,barparams,barpolycoeff = galacticbar

    tstrippy.integrator.setstaticgalaxy(MWname,MWparams)
    tstrippy.integrator.setbackwardorbit()
    tstrippy.integrator.setintegrationparameters(currenttime.value,dt.value,Ntimestep)
    tstrippy.integrator.setinitialkinematics(x0,y0,z0,vx0,vy0,vz0)
    tstrippy.integrator.initgalacticbar(barname,barparams,barpolycoeff)
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