#!/usr/bin/env python
"""
Gravity-coupled orbit integration test.

Demonstrates:
1. Setting up gravity with a single Plummer sphere component
2. Initializing simulator with particles
3. Running leapfrog integration with gravitational forces
4. Checking energy conservation as a sanity check
"""

import numpy as np
import tstrippy
import h5py

def run_gravity_integration():
    """Run a simple N-body integration with a Plummer sphere potential."""
    
    np.random.seed(42)
    
    # =========================================================================
    # STEP 1: Setup gravity module
    # =========================================================================
    print("=" * 70)
    print("STEP 1: Initialize Gravity Module (Plummer Sphere)")
    print("=" * 70)
    
    g = tstrippy.gravity
    g.clear()
    g.set_gravitational_constant(1.0)  # Use natural/physical units (not Galactic)
    
    # Add a Plummer sphere with M=1.0, a=1.0
    # Force: F = -GM*r / (r^2 + a^2)^(3/2)
    M_plummer = 1.0
    a_plummer = 1.0
    g.add_component("plummer", [M_plummer, a_plummer])
    g.finalize()
    
    print(f"✓ Plummer sphere added: M={M_plummer}, a={a_plummer}")
    print(f"  GRAVITY_FINALIZED = {bool(g.gravity_finalized)}")
    print(f"  GRAVITY_NCOMP = {int(g.gravity_ncomp)}")
    
    # =========================================================================
    # STEP 2: Setup simulator with test particles
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 2: Initialize Simulator")
    print("=" * 70)
    
    sim = tstrippy.simulator
    sim.clear()
    
    # Create particles in a circular orbit around the Plummer sphere
    # For circular orbit in a Plummer potential: v_circ = sqrt(GM / (r^2 + a^2)^(1/2))
    n_particles = 20
    r_orbit = 5.0  # orbital radius
    
    # Place particles on a circle in the xy-plane
    theta = np.linspace(0, 2*np.pi, n_particles, endpoint=False)
    x = r_orbit * np.cos(theta)
    y = r_orbit * np.sin(theta)
    z = np.zeros(n_particles)
    
    # Circular velocity: v_c^2 = GM * r / (r^2 + a^2)^(3/2)
    denom = (r_orbit**2 + a_plummer**2)**(1.5)
    v_circ = np.sqrt(M_plummer * r_orbit / denom)
    
    # Velocity perpendicular to radius (counterclockwise)
    vx = -v_circ * np.sin(theta)
    vy =  v_circ * np.cos(theta)
    vz = np.zeros(n_particles)
    
    print(f"✓ Created {n_particles} test particles in circular orbit")
    print(f"  Orbital radius: {r_orbit}")
    print(f"  Circular velocity: {v_circ:.6f}")
    print(f"  Initial kinetic energy per particle: {0.5 * v_circ**2:.6f}")
    
    # Set initial conditions
    sim.setinitialconditions(x, y, z, vx, vy, vz)
    
    # Setup integration: leapfrog with dt=0.01, 500 steps = 5 time units
    dt = 0.01
    nsteps = 500
    sim.setscheme("leapfrog", [0.0, dt, nsteps])
    
    # Enable snapshot writing every 50 steps
    sim.initwritesnapshots(50, "./snapshots", "snapshots")
    
    # Finalize simulator
    sim.finalize()
    
    print(f"\n✓ Simulator configured:")
    print(f"  Scheme: leapfrog")
    print(f"  dt: {dt}")
    print(f"  nsteps: {nsteps}")
    print(f"  Total integration time: {dt * nsteps}")
    
    # =========================================================================
    # STEP 3: Run integration
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 3: Run Integration with Gravity")
    print("=" * 70)
    
    sim.run()
    
    # =========================================================================
    # STEP 4: Analyze results
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 4: Analysis & Energy Conservation Check")
    print("=" * 70)
    
    # Get final state
    x_final = np.asarray(sim.x)
    y_final = np.asarray(sim.y)
    z_final = np.asarray(sim.z)
    vx_final = np.asarray(sim.vx)
    vy_final = np.asarray(sim.vy)
    vz_final = np.asarray(sim.vz)
    
    # Compute final kinetic energy
    v2_final = vx_final**2 + vy_final**2 + vz_final**2
    KE_final = 0.5 * v2_final
    
    # Compute potential energy using gravity module
    # Plummer potential: Phi = -GM / sqrt(r^2 + a^2)
    r2_final = x_final**2 + y_final**2 + z_final**2
    phi_final = -M_plummer / np.sqrt(r2_final + a_plummer**2)
    
    # Total energy (per particle)
    E_final = KE_final + phi_final
    
    # Initial kinetic energy
    KE_initial = 0.5 * (vx**2 + vy**2 + vz**2)
    r2_initial = x**2 + y**2 + z**2
    phi_initial = -M_plummer / np.sqrt(r2_initial + a_plummer**2)
    E_initial = KE_initial + phi_initial
    
    # Energy change
    dE = E_final - E_initial
    dE_rel = np.abs(dE) / np.abs(E_initial)
    
    print(f"\nInitial orbital parameters:")
    print(f"  KE (per particle): {KE_initial[0]:.8f}")
    print(f"  PE (per particle): {phi_initial[0]:.8f}")
    print(f"  E  (per particle): {E_initial[0]:.8f}")
    
    print(f"\nFinal orbital parameters:")
    print(f"  KE (per particle): {KE_final[0]:.8f}")
    print(f"  PE (per particle): {phi_final[0]:.8f}")
    print(f"  E  (per particle): {E_final[0]:.8f}")
    
    print(f"\nEnergy conservation:")
    print(f"  ΔE (per particle): {dE[0]:.2e}")
    print(f"  ΔE / |E| (rel):   {dE_rel[0]:.2e}")
    
    max_dE_rel = np.max(dE_rel)
    print(f"  Max ΔE / |E| (all particles): {max_dE_rel:.2e}")
    
    if max_dE_rel < 1e-4:
        print("  ✓ PASS: Energy conserved to better than 0.01%")
    elif max_dE_rel < 1e-3:
        print("  ⚠ WARN: Energy conserved to 0.1% (acceptable for leapfrog)")
    else:
        print("  ✗ FAIL: Energy conservation poor (check dt or integrator)")
    
    # Check orbital radius change
    r_final = np.sqrt(x_final**2 + y_final**2 + z_final**2)
    r_initial = np.sqrt(x**2 + y**2 + z**2)
    dr = r_final - r_initial
    
    print(f"\nOrbital radius change:")
    print(f"  Initial r: {r_initial[0]:.6f}")
    print(f"  Final r:   {r_final[0]:.6f}")
    print(f"  Δr:        {dr[0]:.2e}")
    print(f"  Δr / r:    {np.abs(dr[0]) / r_initial[0] * 100:.4f}%")
    
    # =========================================================================
    # STEP 5: Timing summary
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 5: Timing Summary")
    print("=" * 70)
    
    print(f"  Total run time:     {sim.timer_run_seconds:.6f} s")
    print(f"  Integration time:   {sim.timer_scheme_seconds:.6f} s ({sim.timer_scheme_seconds/sim.timer_run_seconds*100:.1f}%)")
    print(f"  I/O time:           {sim.timer_write_snapshots_seconds:.6f} s ({sim.timer_write_snapshots_seconds/sim.timer_run_seconds*100:.1f}%)")
    print(f"  Finalize time:      {sim.timer_finalize_seconds:.6f} s")
    
    # =========================================================================
    # STEP 6: Write HDF5 output
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 6: Writing HDF5 Output")
    print("=" * 70)
    
    output_file = "./gravity_integration_example.h5"
    tstrippy.io.write_simulation_hdf5(sim, output_file)
    
    with h5py.File(output_file, "r") as f:
        print(f"\n✓ Wrote {output_file}")
        print(f"  Top-level groups: {list(f.keys())}")
        print(f"  Snapshots saved: {int(f['snapshots'].attrs['nsaved'])}")
        print(f"  Particles tracked: {int(f.attrs['nparticles'])}")
    
    return sim, g, E_initial, E_final


if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("GRAVITY-COUPLED ORBIT INTEGRATION TEST")
    print("=" * 70)
    
    sim, g, E_init, E_final = run_gravity_integration()
    
    print("\n" + "=" * 70)
    print("✓ TEST COMPLETE")
    print("=" * 70)
