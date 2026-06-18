# Bessel Backend: Profile-Agnosticism Refactor

## Summary
Successfully removed the hidden parameter-layout assumptions from the bessel backend. The module is now decoupled from density profile parameter layouts when choosing table-domain scales.

Important scope note:
- This refactor addresses parameter-layout agnosticism only.
- It does **not** prove that the current Bessel table builder is a fully profile-agnostic axisymmetric Poisson solver in the mathematical sense.
- The current implementation still makes solver-level assumptions through its table domain, vertical symmetry handling, and use of a surface-density-style Hankel kernel.

## Changes Made

### 1. **Removed Hardcoded Parameter Extraction**
**Before:**
```fortran
IF (SIZE(params) < 3) THEN
    WRITE(*,'(A)') "WARNING: build_bessel_table_from_density: expected at least 3 params"
    RETURN
END IF

r_scale = MAX(params(2), 1.0D-8)
z_scale = MAX(params(3), 1.0D-8)
```

**After:**
```fortran
! Use public module-level domain scales (profile-agnostic)
r_scale = MAX(BESSEL_R_SCALE, 1.0D-8)
z_scale = MAX(BESSEL_Z_SCALE, 1.0D-8)
```

**Impact:** The bessel module no longer assumes anything about the params array structure. It treats all parameters as profile-intrinsic quantities and derives nothing from them.

---

### 2. **Exposed All Tunable Settings (Like Spherical Harmonics)**

**Module-level PUBLIC parameters (defaults):**
```fortran
REAL*8, PARAMETER, PUBLIC :: BESSEL_G_DEFAULT = 4.30091727D-6
INTEGER, PARAMETER, PUBLIC :: BESSEL_TABLE_NR_DEFAULT = 128
INTEGER, PARAMETER, PUBLIC :: BESSEL_TABLE_NZ_DEFAULT = 128
INTEGER, PARAMETER, PUBLIC :: NK_BUILD_DEFAULT = 256
REAL*8, PARAMETER, PUBLIC :: BESSEL_R_SCALE_DEFAULT = 1.0D0
REAL*8, PARAMETER, PUBLIC :: BESSEL_Z_SCALE_DEFAULT = 1.0D0
```

**Module-level PUBLIC variables (runtime configurable):**
```fortran
REAL*8, PUBLIC :: BESSEL_G = BESSEL_G_DEFAULT
INTEGER, PUBLIC :: BESSEL_TABLE_NR = BESSEL_TABLE_NR_DEFAULT
INTEGER, PUBLIC :: BESSEL_TABLE_NZ = BESSEL_TABLE_NZ_DEFAULT
INTEGER, PUBLIC :: NK_BUILD = NK_BUILD_DEFAULT
REAL*8, PUBLIC :: BESSEL_R_SCALE = BESSEL_R_SCALE_DEFAULT
REAL*8, PUBLIC :: BESSEL_Z_SCALE = BESSEL_Z_SCALE_DEFAULT
```

**Usage pattern (matches spherical harmonics):**
```python
import gravitymini

# Modify defaults before initialization
gravitymini.gravity.BESSEL_TABLE_NR = 256
gravitymini.gravity.BESSEL_TABLE_NZ = 256
gravitymini.gravity.NK_BUILD = 512
gravitymini.gravity.BESSEL_R_SCALE = 2.0
gravitymini.gravity.BESSEL_Z_SCALE = 3.0

# Then call gravity setup as normal
```

---

### 3. **Updated bessel_clear() to Reset All Settings**

Now resets table resolution, quadrature nodes, and domain scales to defaults:
```fortran
BESSEL_G = BESSEL_G_DEFAULT
BESSEL_TABLE_NR = BESSEL_TABLE_NR_DEFAULT
BESSEL_TABLE_NZ = BESSEL_TABLE_NZ_DEFAULT
NK_BUILD = NK_BUILD_DEFAULT
BESSEL_R_SCALE = BESSEL_R_SCALE_DEFAULT
BESSEL_Z_SCALE = BESSEL_Z_SCALE_DEFAULT
```

---

## Architectural Benefits

### ✅ **True Profile-Agnosticism**
- No assumptions about params array layout
- Profile params only encode profile-intrinsic quantities
- Table domain controlled independently via module state

### Current limitation
- The backend is profile-agnostic with respect to parameter encoding, but physics validation is still pending.
- In particular, the current build path integrates the density over `z` into `Sigma(R)` before constructing the Bessel kernel, so thick-disk correctness must be checked explicitly against a 3D reference.

### ✅ **User Control**
- All numerical configuration exposed as public variables
- Can override before initialization (same pattern as spherical harmonics)
- Sensible defaults (r_scale=1.0, z_scale=1.0) for generic domains

### ✅ **Cleaner Separation of Concerns**
- **Profile code**: Provides density via callback with profile params
- **Bessel method**: Uses generic table domain scales, not profile-specific data

---

## Table Domain Interpretation

The bessel table spans:
- **Radial (log scale):** $\log(10^{-3} \times r_{\text{scale}})$ to $\log(10^2 \times r_{\text{scale}})$
- **Vertical (linear):** $0$ to $10^2 \times z_{\text{scale}}$

### Default Domain (r_scale=1.0, z_scale=1.0):
- R: $10^{-3}$ to $10^2$
- Z: $0$ to $10^2$
- Grid: 128 × 128 (bicubic Hermite, 256-node Gauss-Legendre quadrature)

### User Override Example:
```python
# For a galaxy with scale-free structure spanning wider range
BESSEL_R_SCALE = 5.0      # R: 5×10⁻³ to 5×10²
BESSEL_Z_SCALE = 2.0      # Z: 0 to 2×10²
BESSEL_TABLE_NR = 256     # Finer radial resolution
BESSEL_TABLE_NZ = 256     # Finer vertical resolution
```

---

## Validation

✅ **Build:** Successful (Fortran + f2py wrapping)
✅ **Tests:** All 18 smoke tests pass
✅ **Backwards Compatibility:** No changes to public API
✅ **Code Review:** Params extraction removed, domain scales are now profile-agnostic
⚠️ **Physics Status:** Table-build numerics are not yet validated for full thick-disk correctness

---

## Next Steps

1. Run explicit physics gates: symmetry, force/potential finite differences, far-field normalization
2. Add a slow thick-disk reference test for `exponentialdisk`
3. Determine whether the current kernel is a thin-disk / vertically-collapsed approximation or a valid 3D solve
4. Only after physics validation, continue user-guide documentation for recommended settings

## Review Outcome

Current best hypothesis from code review:

- The refactor succeeded at removing hidden `params(2:3)` extraction.
- The likely source of the remaining physics issue is not the profile-agnostic refactor itself.
- The likely issue is in the mathematical form of the table build: it first computes `Sigma(R)` and then reconstructs `Phi(R,z)` with `exp(-k z)`, which is consistent with a vertically-collapsed treatment and therefore must be validated carefully for a genuinely thick density law.

First concrete failure observed after review:

- A new far-field normalization test fails strongly.
- For `Sigma0=1`, `hR=4`, `hZ=0.8` at `R=80`, the backend potential is orders of magnitude larger than the expected monopole limit from `M_tot = 2*pi*Sigma0*hR^2`.
- This strongly suggests that the remaining issue is mathematical / normalization related inside the table build, not API dispatch wiring.

---

## Files Modified

- [besselbfe.f90](besselbfe.f90): Module header and build routine
