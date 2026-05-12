#!/usr/bin/env python
"""
Quick test: Verify exponential_disk_bessel model is registered and can be added.
"""
import sys
sys.path.insert(0, '/Users/sferrone/repos/tstrippy/TEMP')

import tstrippy
import numpy as np

print("=" * 70)
print("Testing exponential_disk_bessel registration and basic setup")
print("=" * 70)

g = tstrippy.gravity
g.clear()

# Try to add exponential_disk_bessel component
# exponentialdisk takes 3 params: Sigma0, scale_r (hR), scale_z (hZ)
params = np.array([1.0, 4.0, 4.0], dtype=np.float64)

try:
    g.add_component("exponentialdisk", params)
    print("✓ exponentialdisk (bessel backend) component added successfully")
except Exception as e:
    print(f"✗ FAILED to add component: {e}")
    raise

try:
    g.finalize()
    print("✓ finalizegravity() succeeded")
except Exception as e:
    print(f"✗ FAILED to finalize: {e}")
    raise

# Try to evaluate at a few points
x = np.array([1.0, 2.0], dtype=float)
y = np.array([0.0, 0.0], dtype=float)
z = np.array([0.0, 0.1], dtype=float)

try:
    ax, ay, az = g.force(x, y, z)
    print(f"✓ Force evaluation succeeded: ax={ax}")
except Exception as e:
    print(f"✗ Force evaluation failed: {e}")
    raise

try:
    phi = g.potential(x, y, z)
    print(f"✓ Potential evaluation succeeded: phi={phi}")
except Exception as e:
    print(f"✗ Potential evaluation failed: {e}")
    raise

print("\n" + "=" * 70)
print("SUCCESS: Bessel backend integrated into gravity dispatch!")
print("=" * 70)
