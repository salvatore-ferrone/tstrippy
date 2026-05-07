import gravitymini
import numpy as np 
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 11,})

r0 = 14
### set the grid points
npoints = 110
x = np.linspace(-2*r0,2*r0,npoints); z = np.linspace(-2*r0,2*r0,npoints)
xq = np.linspace(-2*r0,2*r0,npoints//5); zq = np.linspace(-2*r0,2*r0,npoints//5)
X,Z = np.meshgrid(x,z,indexing="xy")
Xq,Zq = np.meshgrid(xq,zq,indexing="xy")
xq=Xq.flatten() ; zq=Zq.flatten() 
coordsQ=(xq,np.zeros_like(xq),zq)
coordsPHI = (X.flatten(),np.zeros_like(X.flatten()),Z.flatten())


grav = gravitymini.gravity
grav.cleargravity()
grav.addgravitycomponent("exponentialdisk", [1e6, 4, 4])
grav.finalizegravity()

ax,_,az=grav.force(xq,np.zeros_like(xq),zq)
phi=grav.potential(X.flatten(),np.zeros_like(X.flatten()),Z.flatten())
phi = np.reshape(phi,X.shape)

fig,axis=plt.subplots(1,1)
axis.pcolormesh(X,Z,phi)
axis.quiver(xq,zq,ax,az)
axis.set(aspect="equal",title="exponentialdisk")
fig.savefig("ok.png",dpi=300)


grav.cleargravity()
grav.addgravitycomponent("exponentialdisk", [1e6, 4, 4])
grav.finalizegravity()
