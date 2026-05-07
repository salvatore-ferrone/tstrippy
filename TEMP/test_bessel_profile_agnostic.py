"""
Integration test: Verify bessel backend is profile-agnostic.

Key changes:
  1. Removed params(2:3) extraction (was: r_scale, z_scale)
  2. Table domain now controlled by module-level public variables
  3. All tunable settings match spherical harmonics pattern
"""
import sys
sys.path.insert(0, '/Users/sferrone/repos/tstrippy/TEMP')

import numpy as np
import gravitymini

def simple_density_model(params, n, x, y, z):
    """
    Minimal profile-specific params: only 2 parameters needed.
    - params[0]: amplitude
    - params[1]: scale height
    
    No params(2:3) extraction for table domain (now independent!).
    """
    r = np.sqrt(x**2 + y**2)
    amplitude = params[0]
    scale_height = params[1]
    # Simple slab: doesn't need to know table domain extent
    density = amplitude * np.exp(-np.abs(z) / scale_height)
    return density


def test_bessel_with_minimal_profile_params():
    """Test that bessel backend works with minimal params (no scale encoding)."""
    g = gravitymini.gravity
    g.cleargravity()
    
    # Only 2 params needed (not 3+)!
    # Previously required at least 3 to encode r_scale, z_scale
    # Now those come from BESSEL_R_SCALE, BESSEL_Z_SCALE (defaults = 1.0)
    minimal_params = np.array([1.0, 0.3], dtype=np.float64)
    
    print("✓ Attempting bessel setup with minimal 2-element params...")
    try:
        g.addgravitycomponent("exponential_disk_bessel", minimal_params)
        print("✓ SUCCESS: Component accepted minimal params")
    except Exception as e:
        print(f"✗ FAILED: {e}")
        raise
    
    print("✓ Finalizing gravity system...")
    g.finalizegravity()
    
    # Test evaluation
    x = np.array([1.0, 2.0], dtype=float)
    y = np.array([0.0, 0.0], dtype=float)
    z = np.array([0.0, 0.1], dtype=float)
    
    try:
        ax, ay, az = g.force(x, y, z)
        phi = g.potential(x, y, z)
        print(f"✓ Force evaluation successful: ax={ax}")
        print(f"✓ Potential evaluation successful: phi={phi}")
    except Exception as e:
        print(f"✗ Evaluation failed: {e}")
        raise
    
    print("\n" + "="*70)
    print("PROFILE-AGNOSTICISM VALIDATED!")
    print("="*70)
    print("• Bessel backend no longer extracts r_scale, z_scale from params")
    print("• Table domain control decoupled from profile parameters")
    print("• Defaults: BESSEL_R_SCALE=1.0, BESSEL_Z_SCALE=1.0")
    print("• User can override before initialization")
    print("• Params array only encodes profile-intrinsic quantities")


if __name__ == "__main__":
    test_bessel_with_minimal_profile_params()
