from temp import temp as temp 
import galpy.df
import galpy.potential
import numpy as np 
import matplotlib.pyplot as plt
import datetime

# pick inital conditions 
W0 = 5
M=1
rt=1
nparticles=10
npoints_ODE = 1000
nstep = 500000
nskip = 1000
dtfactor = 1e-5

# generate particles thanks to bovy
myking=galpy.df.kingdf(W0=W0, M=M, rt=rt)
R,vR,vT,z,vz,phi=myking.sample(n=nparticles,return_orbit=False)

# put them in cartesian coordinates
xp,yp,zp=R*np.cos(phi),R*np.sin(phi),z
vxp,vyp,vzp=vR*np.cos(phi)-vT*np.sin(phi),vR*np.sin(phi)+vT*np.cos(phi),vz

# compute the dynamic time
v = np.sqrt(vxp**2 + vyp**2 + vzp**2)
r = np.sqrt(xp**2 + yp**2 + zp**2)
tdyn = r/v
dt=np.median(tdyn)*dtfactor

# load the solution from galpy
r,W,dwdr=myking._scalefree_kdf._r,myking._scalefree_kdf._W,myking._scalefree_kdf._dWdr
temp.load_unscaled_king(r,W,dwdr)

# integrate hamilton's equations
starttime = datetime.datetime.now()
temp.leapfrog(dt,nstep,xp,yp,zp,vxp,vyp,vzp,nskip)
endtime = datetime.datetime.now()
print("Integration time: ",endtime-starttime)
time_per_particle_per_step = (endtime-starttime).total_seconds()/(nparticles*nstep)
print("Time per particle per step: ",time_per_particle_per_step)
xt0 = temp.xt
yt0 = temp.yt
zt0 = temp.zt
vxt0 = temp.vxt
vyt0 = temp.vyt
vzt0 = temp.vzt
ntimesamples = xt0.shape[1]
t = nskip*np.arange(ntimesamples)*dt

# store some key quantities
tidal_radius = temp.tidal_radius_
mass = temp.scalefree_mass_
phi0 = -mass/tidal_radius

def mypotential(r):
    if r > tidal_radius:
        return -mass/r
    else:
        return - np.interp(r,myking._scalefree_kdf._r,myking._scalefree_kdf._W) + phi0
    
# compute the energy
rt0 = np.sqrt(xt0**2 + yt0**2 + zt0**2)
pot = np.zeros_like(rt0)
Kinetic = np.zeros_like(pot)
Energy = np.zeros_like(pot)
deltaE = np.zeros_like(pot)
for i in range(nparticles):
    Kinetic[i] = 0.5*(vxt0[i]**2 + vyt0[i]**2 + vzt0[i]**2)
    for j in range(ntimesamples):
        pot[i][j]=mypotential(rt0[i][j])
    Energy[i] = Kinetic[i] + pot[i]
deltaE = np.abs(np.sqrt((Energy - Energy[:,0][:,None])**2)/Energy[:,0][:,None])

fig,axis=plt.subplots(1,2,figsize=(10,5))
for i in range(nparticles):
    axis[0].plot(t,deltaE[i])
    axis[1].plot(xt0[i],yt0[i],alpha=0.5)
axis[1].set_aspect('equal')
axis[1].set_xlabel('x')
axis[1].set_ylabel('y')
axis[0].set_xlabel('Time')
axis[0].set_ylabel('Relative Energy')
axis[0].set_yscale('log')
fig.savefig("energyerror.png")