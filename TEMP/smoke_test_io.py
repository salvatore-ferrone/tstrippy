import tstrippy
import numpy as np

# Check new simulator public flags exist before run
print('simulator.run_success (pre-run):', tstrippy.simulator.run_success)
print('simulator.did_write_snapshots (pre-run):', tstrippy.simulator.did_write_snapshots)
print('simulator.did_write_orbits (pre-run):', tstrippy.simulator.did_write_orbits)

# Quick run with a gravity component
x = np.random.rand(5)
y = x.copy(); z = x.copy()
vx = x.copy(); vy = x.copy(); vz = x.copy()

tstrippy.simulator.clear()
tstrippy.gravity.clear()
tstrippy.simulator.setinitialconditions(x, y, z, vx, vy, vz)
tstrippy.simulator.setscheme('leapfrog', [0.0, 0.01, 10])
tstrippy.simulator.finalize()
tstrippy.simulator.run()

print('After run - run_success:', tstrippy.simulator.run_success)
print('After run - did_write_snapshots:', tstrippy.simulator.did_write_snapshots)
print('After run - did_write_orbits:', tstrippy.simulator.did_write_orbits)

# Check gravity accessors
ncomp = tstrippy.gravity.gravity_ncomp
print('gravity ncomp:', ncomp)
for i in range(1, ncomp + 1):
    name = tstrippy.gravity.getcomponentmodelname(i).decode('ascii').strip()
    npar = int(tstrippy.gravity.getcomponentnparams(i))
    params = tstrippy.gravity.gravity_params[:npar, i - 1]
    print(f'  component {i}: model={name!r}  nparams={npar}  params={params}')

# Check io exports
print('io.write_simulation_hdf5:', tstrippy.io.write_simulation_hdf5)
print('io.read_snapshot_binary:', tstrippy.io.read_snapshot_binary)
print('io.read_orbit_binary:', tstrippy.io.read_orbit_binary)

# Write a test HDF5 file
import tempfile, os
with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as tmp:
    tmpname = tmp.name

tstrippy.io.write_simulation_hdf5(tstrippy.simulator, tstrippy.gravity, tmpname)
import h5py
with h5py.File(tmpname, 'r') as f:
    print('HDF5 groups:', list(f.keys()))
    print('snapshots/state shape:', f['snapshots/state'].shape)
    print('snapshots/time:', f['snapshots/time'][:])
    print('gravity attrs:', dict(f['gravity'].attrs))
    print('gravity/component_01 model_name:', f['gravity/component_01'].attrs['model_name'])
    print('gravity/component_01 params:', f['gravity/component_01/params'][:])

os.unlink(tmpname)
print('SMOKE TEST PASSED')
