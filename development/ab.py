import tstrippy
import galpy.df
import galpy.potential
import numpy as np
import datetime
import matplotlib.pyplot as plt
# let us try integrating a plummer sphere to see the energy and how long it takes 

# pick inital conditions
M=1
b=1
nparticles=100
dtimefactor = 1e-5
nsteps=5e5

nskip = 1000
# generate particles thanks to bovy
pot = galpy.potential.PlummerPotential(amp=M,b=b)
myplummer=galpy.df.isotropicPlummerdf(pot=pot)
R,vR,vT,z,vz,phi=myplummer.sample(n=nparticles,return_orbit=False)
xp,yp,zp=R*np.cos(phi),R*np.sin(phi),z
vxp,vyp,vzp=vR*np.cos(phi)-vT*np.sin(phi),vR*np.sin(phi)+vT*np.cos(phi),vz

# compute the dynamic time
v = np.sqrt(vxp**2 + vyp**2 + vzp**2)
r = np.sqrt(xp**2 + yp**2 + zp**2)
tdyn = r/v
dt=np.median(tdyn)*dtimefactor


# for the static galaxy
staticgalaxy = [1,M,b] # G,M,b
potentialname = "plummer"
# for the integration parameters 
setintegrationparameters = [0.0,dt,nsteps]
# for the initial kinematics
initialkinematics = [xp,yp,zp,vxp,vyp,vzp]
tstrippy.integrator.setstaticgalaxy(potentialname,staticgalaxy)
tstrippy.integrator.setintegrationparameters(*setintegrationparameters)
tstrippy.integrator.setinitialkinematics(*initialkinematics)
# integrate hamilton's equations
starttime   =   datetime.datetime.now()
xt,yt,zt,vxt,vyt,vzt=tstrippy.integrator.leapfrogintime(nsteps,nparticles)
endtime     =   datetime.datetime.now()

print("Integration time: ",endtime-starttime)
time_per_particle_per_step = (endtime-starttime).total_seconds()/(nparticles*nsteps)
print("Time per particle per step: ",time_per_particle_per_step)

# down size the arrays 
xt=xt[:,::nskip]
yt=yt[:,::nskip]
zt=zt[:,::nskip]
vxt=vxt[:,::nskip]
vyt=vyt[:,::nskip]
vzt=vzt[:,::nskip]

rt = np.sqrt(xt**2 + yt**2 + zt**2)
vt = np.sqrt(vxt**2 + vyt**2 + vzt**2)
tdyn = rt/vt
t = np.arange(nsteps+1)*dt
t = t[::nskip]

pot = -M/np.sqrt(rt**2 + b**2)
print(pot.shape)
Energy = 0.5*vt**2 + pot
E0 = Energy[:,0]
deltaE = np.zeros_like(Energy)
for i in range(nparticles):
    deltaE[i] = np.abs( np.sqrt((Energy[i] - E0[i])**2)/E0[i])



print("shape deltaE", deltaE.shape)
print("shape t", t.shape)

fig,axis=plt.subplots(1,2,figsize=(10,5))
for i in range(nparticles):
    axis[0].plot(t,deltaE[i])
    axis[1].plot(xt[i],yt[i],alpha=0.5)
axis[1].set_aspect('equal')
axis[1].set_xlabel('x')
axis[1].set_ylabel('y')
axis[0].set_xlabel('Time')
axis[0].set_ylabel('Relative Energy')
axis[0].set_yscale('log')
text="computation time per particle per step: "+str(time_per_particle_per_step)
fig.text(0.5,0.01,text,ha='center',va='center',transform=axis[0].transAxes)
fig.tight_layout()  

fig.suptitle("Plummer sphere energy conservation")

fig.savefig("plummer_energy_conservation.png")