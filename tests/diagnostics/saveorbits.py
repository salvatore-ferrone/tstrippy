import tstrippy
import numpy as np 
import matplotlib.pyplot as plt 



directory_temporary_snapshots = "./snapshots" # where the binary files will be saved
name_temporary_snapshots = "snapshot"
directory_out_file = "./"
name_final = "streamoutput"

npoints = int(1e2)
x=np.random.rand(npoints);  y=np.random.rand(npoints);  z=np.random.rand(npoints)
vx=np.random.rand(npoints);vy=np.random.rand(npoints); vz=np.random.rand(npoints);

tstrippy.simulator.clear()
tstrippy.simulator.setinitialconditions(x,y,z,vx,vy,vz)
tstrippy.simulator.setscheme("leapfrog",[0.0,1e-2,int(1e2)])
tstrippy.simulator.initwritesnapshots(2,"./snapshots", "snapshots")
tstrippy.simulator.finalize()
tstrippy.simulator.run()
# tstrippy.io.writesnapshot(tstrippy.simulator)