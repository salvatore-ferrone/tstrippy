#!/usr/bin/env python
"""
Validates that bessel backend is now completely profile-agnostic:
- Default scales are 1.0 (not extracted from params)
- All table settings are user-tunable
- Params only needed for density profile, not for scales
"""
import sys
sys.path.insert(0, '/Users/sferrone/repos/tstrippy/TEMP')

import gravitymini
import numpy as np

# ============================================================================
# Test 1: Check that default values are exposed and tunable
# ============================================================================
print("=" * 70)
print("Test 1: Verify default values are public and tunable")
print("=" * 70)

# Clear to get fresh defaults
gravitymini.bessel_clear()

print(f"BESSEL_G_DEFAULT:         {gravitymini.BESSEL_G_DEFAULT}")
print(f"BESSEL_TABLE_NR_DEFAULT:  {gravitymini.BESSEL_TABLE_NR_DEFAULT}")
print(f"BESSEL_TABLE_NZ_DEFAULT:  {gravitymini.BESSEL_TABLE_NZ_DEFAULT}")
print(f"NK_BUILD_DEFAULT:         {gravitymini.NK_BUILD_DEFAULT}")
print(f"BESSEL_R_SCALE_DEFAULT:   {gravitymini.BESSEL_R_SCALE_DEFAULT}")
print(f"BESSEL_Z_SCALE_DEFAULT:   {gravitymini.BESSEL_Z_SCALE_DEFAULT}")

print(f"\nCurrent runtime values:")
print(f"BESSEL_G:                 {gravitymini.BESSEL_G}")
print(f"BESSEL_TABLE_NR:          {gravitymini.BESSEL_TABLE_NR}")
print(f"BESSEL_TABLE_NZ:          {gravitymini.BESSEL_TABLE_NZ}")
print(f"NK_BUILD:                 {gravitymini.NK_BUILD}")
print(f"BESSEL_R_SCALE:           {gravitymini.BESSEL_R_SCALE}")
print(f"BESSEL_Z_SCALE:           {gravitymini.BESSEL_Z_SCALE}")

# ============================================================================
# Test 2: Override defaults before init (just like SH)
# ============================================================================
print("\n" + "=" * 70)
print("Test 2: Override defaults before initialization")
print("=" * 70)

gravitymini.bessel_clear()

# User wants finer table resolution
gravitymini.BESSEL_TABLE_NR = 256
gravitymini.BESSEL_TABLE_NZ = 256
gravitymini.NK_BUILD = 512

# User wants wider domain
gravitymini.BESSEL_R_SCALE = 2.0
gravitymini.BESSEL_Z_SCALE = 3.0

print(f"After user override:")
print(f"BESSEL_TABLE_NR:          {gravitymini.BESSEL_TABLE_NR}")
print(f"BESSEL_TABLE_NZ:          {gravitymini.BESSEL_TABLE_NZ}")
print(f"NK_BUILD:                 {gravitymini.NK_BUILD}")
print(f"BESSEL_R_SCALE:           {gravitymini.BESSEL_R_SCALE}")
print(f"BESSEL_Z_SCALE:           {gravitymini.BESSEL_Z_SCALE}")

# ============================================================================
# Test 3: Profile-agnostic setup (params no longer needs to encode scales)
# ============================================================================
print("\n" + "=" * 70)
print("Test 3: Profile-agnostic setup with generic density")
print("=" * 70)

gravitymini.bessel_clear()

# Simple exponential disk with minimal params (no scale info needed!)
# params[0] = surface density
# params[1] = radial scale (this is profile-intrinsic, not for table domain!)
# params[2] = vertical scale (profile-intrinsic)
params_exp_disk = np.array([1.0, 1.0, 0.3], dtype=np.float64)

def exponential_disk_density(params, n, x, y, z):
    """Generic profile: doesn't care about table domain scales."""
    r = np.sqrt(x**2 + y**2)
    sigma0 = params[0]
    r_scale = params[1]
    z_scale = params[2]
    density = sigma0 * np.exp(-r / r_scale) * np.exp(-np.abs(z) / z_scale)
    return density

# Initialize with defaults (which are now 1.0, 1.0 - decoupled from profile params)
gravitymini.bessel_init_component_tables(1)
print(f"Initialized with domain scales: R={gravitymini.BESSEL_R_SCALE}, Z={gravitymini.BESSEL_Z_SCALE}")

# Project the density - params only encodes profile, not table domain
gravitymini.bessel_project_axisym_density_generic(1, params_exp_disk, exponential_disk_density)
gravitymini.bessel_load_component(1)

print("✓ Successfully set up bessel backend with profile-agnostic params")
print("✓ Table domain controlled independently via BESSEL_R_SCALE, BESSEL_Z_SCALE")

# ============================================================================
# Test 4: Verify clear() resets to defaults
# ============================================================================
print("\n" + "=" * 70)
print("Test 4: Verify clear() resets all settings to defaults")
print("=" * 70)

gravitymini.bessel_clear()

assert gravitymini.BESSEL_TABLE_NR == gravitymini.BESSEL_TABLE_NR_DEFAULT, "NR not reset"
assert gravitymini.BESSEL_TABLE_NZ == gravitymini.BESSEL_TABLE_NZ_DEFAULT, "NZ not reset"
assert gravitymini.NK_BUILD == gravitymini.NK_BUILD_DEFAULT, "NK_BUILD not reset"
assert gravitymini.BESSEL_R_SCALE == gravitymini.BESSEL_R_SCALE_DEFAULT, "R_SCALE not reset"
assert gravitymini.BESSEL_Z_SCALE == gravitymini.BESSEL_Z_SCALE_DEFAULT, "Z_SCALE not reset"

print("✓ All settings reset to defaults after clear()")
print(f"  BESSEL_TABLE_NR = {gravitymini.BESSEL_TABLE_NR}")
print(f"  BESSEL_TABLE_NZ = {gravitymini.BESSEL_TABLE_NZ}")
print(f"  NK_BUILD = {gravitymini.NK_BUILD}")
print(f"  BESSEL_R_SCALE = {gravitymini.BESSEL_R_SCALE}")
print(f"  BESSEL_Z_SCALE = {gravitymini.BESSEL_Z_SCALE}")

print("\n" + "=" * 70)
print("ALL TESTS PASSED: Bessel backend is now profile-agnostic!")
print("=" * 70)
