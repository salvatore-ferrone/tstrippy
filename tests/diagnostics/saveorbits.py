import tstrippy
import numpy as np 
import h5py



directory_temporary_snapshots = "./snapshots" # where the binary files will be saved
name_temporary_snapshots = "snapshots"
output_file = "./streamoutput_example.h5"

npoints = int(1e2)
x=np.random.rand(npoints);  y=np.random.rand(npoints);  z=np.random.rand(npoints)
vx=np.random.rand(npoints);vy=np.random.rand(npoints); vz=np.random.rand(npoints);

tstrippy.simulator.clear()
tstrippy.simulator.setinitialconditions(x,y,z,vx,vy,vz)
tstrippy.simulator.setscheme("leapfrog",[0.0,1e-2,int(1e2)])
tstrippy.simulator.initwritesnapshots(2,"./snapshots", "snapshots")
tstrippy.simulator.finalize()
tstrippy.simulator.run()

tstrippy.io.write_simulation_hdf5(tstrippy.simulator, output_file)

with h5py.File(output_file, "r") as f:
	print("Wrote:", output_file)
	print("Top-level groups:", list(f.keys()))
	print("config groups:", list(f["config"].keys()))
	method = f["config/scheme/method"][()].decode("ascii").strip()
	print("scheme method:", method)
	print("scheme parameters:", f["config/scheme/parameters"][:])
	print("backwardorbit:", bool(f["config"].attrs["backwardorbit"]))
	
	# Timing info
	print("\nTiming (seconds):")
	print("  run:", f["meta/timing/run_seconds"][()])
	print("  finalize:", f["meta/timing/finalize_seconds"][()])
	print("  scheme:", f["meta/timing/scheme_seconds"][()])
	print("  write_snapshots:", f["meta/timing/write_snapshots_seconds"][()])
	print("  write_orbits:", f["meta/timing/write_orbits_seconds"][()])
	
	# Machine info
	print("\nMachine info:")
	print("  architecture:", f["meta/machine/architecture"][()].decode("ascii"))
	print("  processor:", f["meta/machine/processor"][()].decode("ascii"))
	print("  system:", f["meta/machine/system"][()].decode("ascii"))
	print("  release:", f["meta/machine/release"][()].decode("ascii"))
	print("  hostname:", f["meta/machine/hostname"][()].decode("ascii"))