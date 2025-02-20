import matplotlib.pyplot as plt
import matplotlib as mpl
from galpy.orbit import Orbit
from galpy.potential import KingPotential
from galpy.df import kingdf
import numpy as np

# generate 100 particles within a king potential and integrate them 
nparticles = 100
W0 = 10
myking = kingdf(W0)
R,vR,vT,z,vz,phi=myking.sample(n=nparticles,return_orbit=False)
# compute the dynamic time 
v = np.sqrt(vR**2 + vT**2 + vz**2)
r = np.sqrt(R**2 + z**2)
tdynam = r/v
dt = 1e-3*tdynam.min()
nstep = 3000
ttotal = nstep*dt
t= np.arange(nstep+1)*dt
o=Orbit(vxvv=[R,vR,vT,z,vz,phi])
mypotential = KingPotential(ro=1.,vo=1.,W0=W0)
o.integrate(t=t, pot=mypotential)

# get the initial energy from each partcile
E0 = np.zeros(nparticles)
deltaE=np.zeros((nstep+1,nparticles))
for i in range(nparticles):
    E0[i] = o.E(t)[i][0]
    deltaE[:,i] = np.sqrt((o.E(t)[i] - E0[i])**2)/E0[i]
# gget the error


# plot the results
fig,axis=plt.subplots(1,1,figsize=(5,6))
for i in range(nparticles):
    axis.plot(t,deltaE[:,i],alpha=0.5)
axis.set_xlabel("Time [Gyr]")
axis.set_ylabel("Error in Energy [km^2/s^2]")
fig.tight_layout()

fig.savefig("galpy_king_orbits_energy_error_python.png")
plt.close(fig)

# now plot the x-y projection of the orbits

fig,axis=plt.subplots(1,1,figsize=(5,6))
for i in range(nparticles):
    axis.plot(o.x(t)[i],o.y(t)[i],alpha=0.5)
axis.set_xlabel("x [unit]")
axis.set_ylabel("y [unit]")
axis.set_aspect('equal')
fig.tight_layout()
fig.savefig("galpy_king_orbits_xy_python.png")
plt.close(fig)