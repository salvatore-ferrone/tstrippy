#!/usr/bin/env python
"""
Validates that bessel backend now uses spherical-harmonics-style init pattern.

Pattern:
  - bessel_set_gravity_constant(g) - set G separately
  - bessel_init(nr, nz, nk_build, r_scale, z_scale) - explicit init with all params
  - bessel_default_init() - init with hardcoded defaults
  - bessel_clear() - reset to defaults
  - bessel_init_component_tables(ncomp) - allocate component storage (after init)
  - bessel_project_axisym_density_generic(icomp, params, density_func) - fill table
"""
import sys
sys.path.insert(0, '/Users/sferrone/repos/tstrippy/TEMP')

import tstrippy
import numpy as np

# Access bessel backend through besselbfe submodule
bessel = tstrippy.besselbfe

print("=" * 70)
print("Test 1: Verify bessel_default_init() exists and works")
print("=" * 70)

g = tstrippy.gravity
g.cleargravity()

# This should initialize bessel backend with hardcoded defaults
try:
    bessel.bessel_default_init()
    print(f"✓ bessel_default_init() succeeded")
    print(f"  BESSEL_TABLE_NR: {bessel.bessel_table_nr}")
    print(f"  BESSEL_TABLE_NZ: {bessel.bessel_table_nz}")
    print(f"  NK_BUILD: {bessel.nk_build}")
    print(f"  BESSEL_R_SCALE: {bessel.bessel_r_scale}")
    print(f"  BESSEL_Z_SCALE: {bessel.bessel_z_scale}")
    print(f"  BESSEL_G: {bessel.bessel_g}")
    print(f"  BESSEL_INITIALIZED: {bessel.bessel_initialized}")
except Exception as e:
    print(f"✗ FAILED: {e}")
    raise

print("\n" + "=" * 70)
print("Test 2: Explicit bessel_init with custom parameters")
print("=" * 70)

bessel.bessel_clear()
print(f"After clear: BESSEL_INITIALIZED = {bessel.bessel_initialized}")

# Initialize with custom parameters (matching SH pattern)
try:
    bessel.bessel_init(256, 256, 512, 2.0, 3.0)
    print(f"✓ bessel_init(256, 256, 512, 2.0, 3.0) succeeded")
    print(f"  BESSEL_TABLE_NR: {bessel.bessel_table_nr} (expected 256)")
    print(f"  BESSEL_TABLE_NZ: {bessel.bessel_table_nz} (expected 256)")
    print(f"  NK_BUILD: {bessel.nk_build} (expected 512)")
    print(f"  BESSEL_R_SCALE: {bessel.bessel_r_scale} (expected 2.0)")
    print(f"  BESSEL_Z_SCALE: {bessel.bessel_z_scale} (expected 3.0)")
    print(f"  BESSEL_INITIALIZED: {bessel.bessel_initialized} (expected True)")
    
    assert bessel.bessel_table_nr == 256
    assert bessel.bessel_table_nz == 256
    assert bessel.nk_build == 512
    assert np.isclose(bessel.bessel_r_scale, 2.0)
    assert np.isclose(bessel.bessel_z_scale, 3.0)
    assert bessel.bessel_initialized == True
except Exception as e:
    print(f"✗ FAILED: {e}")
    raise

print("\n" + "=" * 70)
print("Test 3: bessel_set_gravity_constant() works independently")
print("=" * 70)

g_test = 1.0
try:
    bessel.bessel_set_gravity_constant(g_test)
    print(f"✓ bessel_set_gravity_constant({g_test}) succeeded")
    print(f"  BESSEL_G: {bessel.bessel_g} (expected {g_test})")
    assert np.isclose(bessel.bessel_g, g_test)
except Exception as e:
    print(f"✗ FAILED: {e}")
    raise

print("\n" + "=" * 70)
print("Test 4: bessel_clear() resets all settings to defaults")
print("=" * 70)

bessel.bessel_clear()
print(f"After clear():")
print(f"  BESSEL_TABLE_NR: {bessel.bessel_table_nr} (expected 128)")
print(f"  BESSEL_TABLE_NZ: {bessel.bessel_table_nz} (expected 128)")
print(f"  NK_BUILD: {bessel.nk_build} (expected 256)")
print(f"  BESSEL_R_SCALE: {bessel.bessel_r_scale} (expected 1.0)")
print(f"  BESSEL_Z_SCALE: {bessel.bessel_z_scale} (expected 1.0)")
print(f"  BESSEL_INITIALIZED: {bessel.bessel_initialized} (expected False)")

assert bessel.bessel_table_nr == 128
assert bessel.bessel_table_nz == 128
assert bessel.nk_build == 256
assert np.isclose(bessel.bessel_r_scale, 1.0)
assert np.isclose(bessel.bessel_z_scale, 1.0)
assert bessel.bessel_initialized == False

print("\n" + "=" * 70)
print("SUCCESS: Bessel backend now uses SH-style init pattern!")
print("=" * 70)
print("Usage pattern:")
print("  1. bessel.bessel_set_gravity_constant(g)")
print("  2. bessel.bessel_init(nr, nz, nk, r_scale, z_scale)")
print("     OR bessel.bessel_default_init()")
print("  3. bessel.bessel_init_component_tables(ncomp)")
print("  4. bessel.bessel_project_axisym_density_generic(...)")
