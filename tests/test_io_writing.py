#!/usr/bin/env python
"""Test snapshot and orbit file writing."""

import os
import tempfile
import numpy as np
import tstrippy

sim = tstrippy.simulator
basecomponent = ['plummer', [1e6, 5e-3]]

def test_snapshot_writing():
    """Test that snapshot files are created during run()."""
    sim.clear()
    
    # Simple 2-particle setup
    x = np.array([0.0, 1.0], dtype=np.float64)
    y = np.array([0.0, 0.0], dtype=np.float64)
    z = np.array([0.0, 0.0], dtype=np.float64)
    vx = np.array([1.0, -1.0], dtype=np.float64)
    vy = np.array([0.0, 0.0], dtype=np.float64)
    vz = np.array([0.0, 0.0], dtype=np.float64)
    
    sim.setinitialconditions(x, y, z, vx, vy, vz)
    # leapfrog params: [gamma, dt, nsteps]
    sim.setscheme("leapfrog", [0.0, 0.01, 10])
    sim.settimestamps(np.linspace(0.0, 10.0, 11))  # 11 timestamps (0,1,2,...,10)
    
    # Setup snapshot writing with temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        sim.initwritesnapshots(nskip=2, directory=tmpdir, basename="snap")
        sim.add_component(*basecomponent)
        # finalize 
        sim.finalize()
        # Run simulation
        sim.run()
        
        # Check that snapshot files were created (nskip=2, nsteps=10, so steps 2,4,6,8,10)
        # With nskip=2 from 10 steps, we write at istep=2,4,6,8,10 = 5 files
        files = sorted([f for f in os.listdir(tmpdir) if f.startswith("snap_")])
        print(f"Snapshot files created: {files}")
        
        # Verify file count
        assert len(files) > 0, "No snapshot files created!"
        
        # Verify each file is readable and has expected structure
        for fname in files:
            fpath = os.path.join(tmpdir, fname)
            assert os.path.getsize(fpath) > 0, f"Snapshot file {fname} is empty!"
            print(f"  {fname}: {os.path.getsize(fpath)} bytes")


def test_orbit_writing():
    """Test that orbit files are created after run()."""
    sim.clear()
    sim.max_ram_mb = 1024.0  # Ensure memory available
    sim.orbit_ram_limit_mb = 1024.0
    
    # Simple 3-particle setup with orbit tracking
    x = np.array([0.0, 1.0, 2.0], dtype=np.float64)
    y = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    z = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    vx = np.array([1.0, 0.5, 0.2], dtype=np.float64)
    vy = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    vz = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    
    sim.setinitialconditions(x, y, z, vx, vy, vz)
    sim.setscheme("leapfrog", [0.0, 0.01, 5])
    sim.add_component(*basecomponent)
    # Track all 3 particles
    sim.trim_orbits(2)  # Skip every 2 timesteps
    
    sim.settimestamps(np.linspace(0.0, 5.0, 6))  # 6 timestamps (0,1,2,3,4,5)
    
    # Setup orbit writing
    with tempfile.TemporaryDirectory() as tmpdir:
        sim.initwriteorbits(nskip=1, directory=tmpdir, basename="orbits")
        
        # Run simulation
        sim.run()
        assert sim.run_success, "the run should be successful"
        # Check that orbit files were created (one per tracked particle)
        files = sorted([f for f in os.listdir(tmpdir) if f.startswith("orbits_particle_")])
        print(f"Orbit files created: {files}")
        
        # We track all 3 particles (trim_orbits set n_particles_orbit to Nparticles)
        # So we expect 3 orbit files
        assert len(files) == 3, f"Expected 3 orbit files, got {len(files)}"
        
        # Verify each file is non-empty
        for fname in files:
            fpath = os.path.join(tmpdir, fname)
            size = os.path.getsize(fpath)
            assert size > 0, f"Orbit file {fname} is empty!"
            print(f"  {fname}: {size} bytes")


def test_snapshot_and_orbit_together():
    """Test simultaneous snapshot and orbit writing."""
    sim.clear()
    
    x = np.array([0.0, 1.0], dtype=np.float64)
    y = np.array([0.0, 0.0], dtype=np.float64)
    z = np.array([0.0, 0.0], dtype=np.float64)
    vx = np.array([1.0, -1.0], dtype=np.float64)
    vy = np.array([0.0, 0.0], dtype=np.float64)
    vz = np.array([0.0, 0.0], dtype=np.float64)
    
    sim.setinitialconditions(x, y, z, vx, vy, vz)
    sim.setscheme("leapfrog", [0.0, 0.01, 4])
    sim.settimestamps(np.linspace(0.0, 4.0, 5))  # 5 timesteps
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Enable both snapshot and orbit writing
        snap_dir = os.path.join(tmpdir, "snapshots")
        orbit_dir = os.path.join(tmpdir, "orbits")
        os.makedirs(snap_dir)
        os.makedirs(orbit_dir)
        
        sim.initwritesnapshots(nskip=1, directory=snap_dir, basename="snapshot")
        sim.initwriteorbits(nskip=1, directory=orbit_dir, basename="orbit")
        sim.add_component(*basecomponent)
        # Run simulation
        sim.run()
        assert sim.run_success, "the run should be successful"
        # Check snapshots
        snap_files = sorted([f for f in os.listdir(snap_dir) if f.startswith("snapshot_")])
        print(f"Snapshots: {snap_files}")
        assert len(snap_files) > 0, "No snapshot files created"
        
        # Check orbits
        orbit_files = sorted([f for f in os.listdir(orbit_dir) if f.startswith("orbit_particle_")])
        print(f"Orbits: {orbit_files}")
        assert len(orbit_files) > 0, "No orbit files created"


if __name__ == "__main__":
    test_snapshot_writing()
    print("✓ Snapshot writing test passed\n")
    
    test_orbit_writing()
    print("✓ Orbit writing test passed\n")
    
    test_snapshot_and_orbit_together()
    print("✓ Combined snapshot+orbit writing test passed\n")
    
    print("All I/O tests passed!")
