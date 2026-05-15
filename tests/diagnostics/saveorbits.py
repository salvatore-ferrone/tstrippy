import tstrippy
import numpy as np 
import matplotlib.pyplot as plt 


npoints = int(1e2)
x=np.random.rand(npoints);y=np.random.rand(npoints);z=np.random.rand(npoints)
vx=np.random.rand(npoints);vy=np.random.rand(npoints);vz=np.random.rand(npoints);

tstrippy.simulator.clear()
tstrippy.simulator.setinitialconditions(x,y,z,vx,vy,vz)
tstrippy.simulator.setscheme("leapfrog",[0.0,1e-2,int(1e2)])
tstrippy.simulator.initwritesnapshots(2,"./snapshots", "snapshots")
tstrippy.simulator.finalize()
tstrippy.simulator.run()
tstrippy.simulator.orbits.shape