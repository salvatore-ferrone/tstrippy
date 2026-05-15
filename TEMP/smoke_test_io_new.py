import tstrippy
import numpy as np
import tempfile
import os

print("=" * 60)
print("SMOKE TEST: I/O Infrastructure (no simulation)")
print("=" * 60)

# 1. Check new simulator public flags exist
print("\n1. Checking simulator public flags...")
assert hasattr(tstrippy.simulator, 'run_success'), "Missing run_success"
assert hasattr(tstrippy.simulator, 'did_write_snapshots'), "Missing did_write_snapshots"
assert hasattr(tstrippy.simulator, 'did_write_orbits'), "Missing did_write_orbits"
print("   ✓ run_success:", tstrippy.simulator.run_success)
print("   ✓ did_write_snapshots:", tstrippy.simulator.did_write_snapshots)
print("   ✓ did_write_orbits:", tstrippy.simulator.did_write_orbits)

# 2. Check gravity state and component model names
print("\n2. Checking gravity state...")
print("   ✓ gravity_ncomp:", tstrippy.gravity.gravity_ncomp)
print("   ✓ gravity_finalized:", tstrippy.gravity.gravity_finalized)
print("   ✓ component_model_names (array):", tstrippy.gravity.component_model_names.shape)

# 3. Check io exports
print("\n3. Checking io module exports...")
assert hasattr(tstrippy.io, 'write_simulation_hdf5'), "Missing write_simulation_hdf5"
assert hasattr(tstrippy.io, 'read_snapshot_binary'), "Missing read_snapshot_binary"
assert hasattr(tstrippy.io, 'read_orbit_binary'), "Missing read_orbit_binary"
print("   ✓ write_simulation_hdf5 callable")
print("   ✓ read_snapshot_binary callable")
print("   ✓ read_orbit_binary callable")

# 4. Test HDF5 file structure with mock data (no run needed)
print("\n4. Testing HDF5 file creation with mock data...")
try:
    import h5py
except ImportError:
    print("   ✗ h5py not installed; skipping file write test")
else:
    with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as tmp:
        tmpname = tmp.name
    
    try:
        # Manually set flags to simulate a completed run
        # (without actually running integration)
        tstrippy.simulator.run_success = 1
        tstrippy.simulator.did_write_snapshots = 0
        tstrippy.simulator.did_write_orbits = 0
        
        # Write HDF5
        tstrippy.io.write_simulation_hdf5(
            tstrippy.simulator,
            tstrippy.gravity,
            tmpname
        )
        
        # Verify structure
        with h5py.File(tmpname, 'r') as f:
            assert 'meta' in f, "Missing 'meta' group"
            assert 'config' in f, "Missing 'config' group"
            assert 'gravity' in f, "Missing 'gravity' group"
            assert 'snapshots' in f, "Missing 'snapshots' group"
            print("   ✓ HDF5 groups created: meta, config, gravity, snapshots")
            
            # Check attributes
            assert 'schema_version' in f['meta'].attrs
            assert 'created_utc' in f['meta'].attrs
            print("   ✓ meta attributes present")
            
            assert 'nparticles' in f['config'].attrs
            print("   ✓ config attributes present")
            
            assert 'G' in f['gravity'].attrs
            assert 'ncomponents' in f['gravity'].attrs
            print("   ✓ gravity attributes present (G, ncomponents)")
            
        # Clean up
        tstrippy.simulator.run_success = 0
        os.unlink(tmpname)
        print("   ✓ HDF5 file successfully created and validated")
        
    except Exception as e:
        print(f"   ✗ Error during HDF5 test: {e}")
        import traceback
        traceback.print_exc()
        if os.path.exists(tmpname):
            os.unlink(tmpname)
        raise

print("\n" + "=" * 60)
print("✓ ALL SMOKE TESTS PASSED")
print("=" * 60)
