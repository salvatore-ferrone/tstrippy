MODULE besselbfe
    USE mathutils, ONLY: gauss_legendre_nodes_weights, bessel_j0_scalar, bessel_j1_scalar, &
                         bicubic_hermite_eval_2d
    IMPLICIT NONE

    ! Private defaults (constants)
    REAL*8, PARAMETER, PRIVATE :: BESSEL_G_DEFAULT = 4.30091727D-6
    INTEGER, PARAMETER, PRIVATE :: BESSEL_TABLE_NR_DEFAULT = 128
    INTEGER, PARAMETER, PRIVATE :: BESSEL_TABLE_NZ_DEFAULT = 128
    INTEGER, PARAMETER, PRIVATE :: NK_BUILD_DEFAULT = 256
    REAL*8, PARAMETER, PRIVATE :: BESSEL_R_SCALE_DEFAULT = 1.0D0
    REAL*8, PARAMETER, PRIVATE :: BESSEL_Z_SCALE_DEFAULT = 1.0D0

    ! Meta field count (internal bookkeeping)
    INTEGER, PARAMETER, PRIVATE :: BESSEL_META_FIELDS = 4

    ! Public state (read-only, set via init routines)
    REAL*8, PUBLIC :: BESSEL_G = BESSEL_G_DEFAULT
    INTEGER, PUBLIC :: BESSEL_TABLE_NR = BESSEL_TABLE_NR_DEFAULT
    INTEGER, PUBLIC :: BESSEL_TABLE_NZ = BESSEL_TABLE_NZ_DEFAULT
    INTEGER, PUBLIC :: NK_BUILD = NK_BUILD_DEFAULT
    REAL*8, PUBLIC :: BESSEL_R_SCALE = BESSEL_R_SCALE_DEFAULT
    REAL*8, PUBLIC :: BESSEL_Z_SCALE = BESSEL_Z_SCALE_DEFAULT

    LOGICAL, PUBLIC :: BESSEL_INITIALIZED = .FALSE.
    INTEGER, PUBLIC :: BESSEL_NCOMP = 0
    INTEGER, PUBLIC :: BESSEL_ACTIVE_COMP = 0

    REAL*8, ALLOCATABLE, PRIVATE :: BESSEL_TABLE_PHI(:, :, :)
    REAL*8, ALLOCATABLE, PRIVATE :: BESSEL_TABLE_DPHI_DR(:, :, :)
    REAL*8, ALLOCATABLE, PRIVATE :: BESSEL_TABLE_DPHI_DZ(:, :, :)
    REAL*8, ALLOCATABLE, PRIVATE :: BESSEL_TABLE_D2PHI_DRDZ(:, :, :)
    REAL*8, ALLOCATABLE, PRIVATE :: BESSEL_TABLE_META(:, :)
    REAL*8, ALLOCATABLE, PRIVATE :: BESSEL_COMPONENT_R_SCALE(:)
    REAL*8, ALLOCATABLE, PRIVATE :: BESSEL_COMPONENT_Z_SCALE(:)
    LOGICAL, ALLOCATABLE, PRIVATE :: BESSEL_COMPONENT_READY(:)

    ABSTRACT INTERFACE
        SUBROUTINE axisymmetric_density_model(params, n, x, y, z, rho)
            REAL*8, INTENT(IN), DIMENSION(:) :: params
            INTEGER, INTENT(IN) :: n
            REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
            REAL*8, INTENT(OUT), DIMENSION(n) :: rho
        END SUBROUTINE axisymmetric_density_model
    END INTERFACE

    PRIVATE :: build_table

CONTAINS

    SUBROUTINE clear()
        IMPLICIT NONE

        BESSEL_G = BESSEL_G_DEFAULT
        BESSEL_TABLE_NR = BESSEL_TABLE_NR_DEFAULT
        BESSEL_TABLE_NZ = BESSEL_TABLE_NZ_DEFAULT
        NK_BUILD = NK_BUILD_DEFAULT
        BESSEL_R_SCALE = BESSEL_R_SCALE_DEFAULT
        BESSEL_Z_SCALE = BESSEL_Z_SCALE_DEFAULT
        BESSEL_INITIALIZED = .FALSE.
        BESSEL_NCOMP = 0
        BESSEL_ACTIVE_COMP = 0

        IF (ALLOCATED(BESSEL_TABLE_PHI)) DEALLOCATE(BESSEL_TABLE_PHI)
        IF (ALLOCATED(BESSEL_TABLE_DPHI_DR)) DEALLOCATE(BESSEL_TABLE_DPHI_DR)
        IF (ALLOCATED(BESSEL_TABLE_DPHI_DZ)) DEALLOCATE(BESSEL_TABLE_DPHI_DZ)
        IF (ALLOCATED(BESSEL_TABLE_D2PHI_DRDZ)) DEALLOCATE(BESSEL_TABLE_D2PHI_DRDZ)
        IF (ALLOCATED(BESSEL_TABLE_META)) DEALLOCATE(BESSEL_TABLE_META)
        IF (ALLOCATED(BESSEL_COMPONENT_R_SCALE)) DEALLOCATE(BESSEL_COMPONENT_R_SCALE)
        IF (ALLOCATED(BESSEL_COMPONENT_Z_SCALE)) DEALLOCATE(BESSEL_COMPONENT_Z_SCALE)
        IF (ALLOCATED(BESSEL_COMPONENT_READY)) DEALLOCATE(BESSEL_COMPONENT_READY)
    END SUBROUTINE clear

    SUBROUTINE set_gravitational_constant(g)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: g

        IF (g <= 0.0D0) THEN
            WRITE(*,'(A)') "WARNING: bessel_set_gravity_constant: g must be positive"
            RETURN
        END IF
        BESSEL_G = g
    END SUBROUTINE set_gravitational_constant

    SUBROUTINE initialize(nr, nz, nk_build_in)
        ! Explicit initialization for table/grid resolution controls.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: nr, nz, nk_build_in

        IF (nr < 2) THEN
            WRITE(*,'(A)') "WARNING: initialize: nr must be >= 2"
            RETURN
        END IF
        IF (nz < 2) THEN
            WRITE(*,'(A)') "WARNING: initialize: nz must be >= 2"
            RETURN
        END IF
        IF (nk_build_in < 4) THEN
            WRITE(*,'(A)') "WARNING: initialize: nk_build must be >= 4"
            RETURN
        END IF
        CALL clear()
        BESSEL_TABLE_NR = nr
        BESSEL_TABLE_NZ = nz
        NK_BUILD = nk_build_in
        BESSEL_R_SCALE = BESSEL_R_SCALE_DEFAULT
        BESSEL_Z_SCALE = BESSEL_Z_SCALE_DEFAULT
        BESSEL_INITIALIZED = .TRUE.
    END SUBROUTINE initialize

    SUBROUTINE default_initialize()
        ! Default initialization using hardcoded parameter defaults
        IMPLICIT NONE
        CALL initialize(BESSEL_TABLE_NR_DEFAULT, BESSEL_TABLE_NZ_DEFAULT, NK_BUILD_DEFAULT)
    END SUBROUTINE default_initialize

    SUBROUTINE allocate_component_tables(ncomp)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ncomp

        IF (ncomp < 1) THEN
            WRITE(*,'(A)') "WARNING: allocate_table: ncomp must be >= 1"
            RETURN
        END IF

        BESSEL_NCOMP = ncomp

        ALLOCATE(BESSEL_TABLE_PHI(BESSEL_TABLE_NR, BESSEL_TABLE_NZ, BESSEL_NCOMP))
        ALLOCATE(BESSEL_TABLE_DPHI_DR(BESSEL_TABLE_NR, BESSEL_TABLE_NZ, BESSEL_NCOMP))
        ALLOCATE(BESSEL_TABLE_DPHI_DZ(BESSEL_TABLE_NR, BESSEL_TABLE_NZ, BESSEL_NCOMP))
        ALLOCATE(BESSEL_TABLE_D2PHI_DRDZ(BESSEL_TABLE_NR, BESSEL_TABLE_NZ, BESSEL_NCOMP))
        ALLOCATE(BESSEL_TABLE_META(BESSEL_META_FIELDS, BESSEL_NCOMP))
        ALLOCATE(BESSEL_COMPONENT_R_SCALE(BESSEL_NCOMP))
        ALLOCATE(BESSEL_COMPONENT_Z_SCALE(BESSEL_NCOMP))
        ALLOCATE(BESSEL_COMPONENT_READY(BESSEL_NCOMP))

        BESSEL_TABLE_PHI = 0.0D0
        BESSEL_TABLE_DPHI_DR = 0.0D0
        BESSEL_TABLE_DPHI_DZ = 0.0D0
        BESSEL_TABLE_D2PHI_DRDZ = 0.0D0
        BESSEL_TABLE_META = 0.0D0
        BESSEL_COMPONENT_R_SCALE = BESSEL_R_SCALE
        BESSEL_COMPONENT_Z_SCALE = BESSEL_Z_SCALE
        BESSEL_COMPONENT_READY = .FALSE.
    END SUBROUTINE allocate_component_tables

    SUBROUTINE set_component_scales(component_index, r_scale, z_scale)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: component_index
        REAL*8, INTENT(IN) :: r_scale, z_scale

        IF (component_index < 1 .OR. component_index > BESSEL_NCOMP) THEN
            WRITE(*,'(A)') "WARNING: set_component_scales: invalid component_index"
            RETURN
        END IF
        IF (r_scale <= 0.0D0) THEN
            WRITE(*,'(A)') "WARNING: set_component_scales: r_scale must be positive"
            RETURN
        END IF
        IF (z_scale <= 0.0D0) THEN
            WRITE(*,'(A)') "WARNING: set_component_scales: z_scale must be positive"
            RETURN
        END IF

        BESSEL_COMPONENT_R_SCALE(component_index) = r_scale
        BESSEL_COMPONENT_Z_SCALE(component_index) = z_scale
    END SUBROUTINE set_component_scales

    SUBROUTINE project_density(component_index, density_model, params)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: component_index
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        PROCEDURE(axisymmetric_density_model) :: density_model
        REAL*8, DIMENSION(1) :: x0, y0, z0, rho0

        IF (.NOT. BESSEL_INITIALIZED) THEN
            WRITE(*,'(A)') "WARNING: project_density: call allocate_table first"
            RETURN
        END IF
        IF (component_index < 1 .OR. component_index > BESSEL_NCOMP) THEN
            WRITE(*,'(A)') "WARNING: project_density: invalid component_index"
            RETURN
        END IF
        IF (SIZE(params) < 1) THEN
            WRITE(*,'(A)') "WARNING: project_density: params must be non-empty"
            RETURN
        END IF

        ! Density callback is part of the generic API contract and used for
        ! quick input sanity checks during projection setup.
        x0(1) = 0.0D0
        y0(1) = 0.0D0
        z0(1) = 0.0D0
        CALL density_model(params, 1, x0, y0, z0, rho0)
        IF (rho0(1) < 0.0D0) THEN
            WRITE(*,'(A)') "WARNING: project_density: density_model returned rho<0 at origin"
        END IF

        CALL build_table(component_index, density_model, params)
        BESSEL_COMPONENT_READY(component_index) = .TRUE.
    END SUBROUTINE project_density

    SUBROUTINE load_component(component_index)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: component_index

        IF (.NOT. BESSEL_INITIALIZED) THEN
            WRITE(*,'(A)') "WARNING: bessel_load_component: backend not initialized"
            RETURN
        END IF
        IF (component_index < 1 .OR. component_index > BESSEL_NCOMP) THEN
            WRITE(*,'(A)') "WARNING: bessel_load_component: invalid component_index"
            RETURN
        END IF
        IF (.NOT. BESSEL_COMPONENT_READY(component_index)) THEN
            WRITE(*,'(A)') "WARNING: bessel_load_component: component table not built"
            RETURN
        END IF

        BESSEL_ACTIVE_COMP = component_index
    END SUBROUTINE load_component

    SUBROUTINE force(n, x, y, z, ax, ay, az)
        ! wrappers
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: ax, ay, az
        REAL*8, DIMENSION(n) :: phi_tmp

        CALL component_force_potential(n, x, y, z, ax, ay, az, phi_tmp)
    END SUBROUTINE force

    SUBROUTINE potential(n, x, y, z, phi)
        ! wrapper 
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: phi
        REAL*8, DIMENSION(n) :: ax_tmp, ay_tmp, az_tmp

        CALL component_force_potential(n, x, y, z, ax_tmp, ay_tmp, az_tmp, phi)
    END SUBROUTINE potential

    SUBROUTINE component_force_potential(n, x, y, z, ax, ay, az, phi)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN),  DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: ax, ay, az, phi
        REAL*8, PARAMETER :: eps = 1.0D-16
        INTEGER :: i, iR, iZ, ic
        REAL*8 :: R, Rq, absz, signz, aR_val
        REAL*8 :: r3, r_sph, mass_est, phi_mono, blend0, blend1, t, w_blend, scale_ref
        REAL*8 :: logR, logR_min, dlogR, dz, logR_max
        REAL*8 :: alphaR, alphaZ, zq
        REAL*8 :: R_lo, R_hi, dR, z_lo
        REAL*8 :: phi_loc, dphi_dR_loc, dphi_dz_abs
        REAL*8 :: f00, f10, f01, f11
        REAL*8 :: gr00, gr10, gr01, gr11
        REAL*8 :: gz00, gz10, gz01, gz11

        ax = 0.0D0
        ay = 0.0D0
        az = 0.0D0
        phi = 0.0D0

        IF (.NOT. BESSEL_INITIALIZED) THEN
            WRITE(*,'(A)') "WARNING: bessel_eval_component: backend not initialized"
            RETURN
        END IF

        ic = BESSEL_ACTIVE_COMP
        IF (ic < 1 .OR. ic > BESSEL_NCOMP) THEN
            WRITE(*,'(A)') "WARNING: bessel_eval_component: no active component loaded"
            RETURN
        END IF
        IF (.NOT. BESSEL_COMPONENT_READY(ic)) THEN
            WRITE(*,'(A)') "WARNING: bessel_eval_component: active component table not built"
            RETURN
        END IF

        logR_min = BESSEL_TABLE_META(1, ic)
        dlogR = BESSEL_TABLE_META(2, ic)
        dz = BESSEL_TABLE_META(3, ic)
        logR_max = logR_min + DBLE(BESSEL_TABLE_NR - 1) * dlogR

        DO i = 1, n
            R = SQRT(x(i)**2 + y(i)**2)
            absz = ABS(z(i))
            IF (z(i) > 0.0D0) THEN
                signz = 1.0D0
            ELSE IF (z(i) < 0.0D0) THEN
                signz = -1.0D0
            ELSE
                signz = 0.0D0
            END IF

            IF (R > eps) THEN
                logR = MIN(MAX(LOG(R), logR_min), logR_max)
                Rq = EXP(logR)
            ELSE
                logR = logR_min
                Rq = EXP(logR_min)
            END IF
            zq = MIN(absz, DBLE(BESSEL_TABLE_NZ - 1) * dz)

            iR = INT((logR - logR_min) / dlogR)
            iR = MAX(0, MIN(BESSEL_TABLE_NR - 2, iR))
            R_lo = EXP(logR_min + DBLE(iR) * dlogR)
            R_hi = EXP(logR_min + DBLE(iR + 1) * dlogR)
            dR = MAX(R_hi - R_lo, eps)
            alphaR = (Rq - R_lo) / dR
            alphaR = MAX(0.0D0, MIN(1.0D0, alphaR))

            iZ = INT(zq / dz)
            iZ = MAX(0, MIN(BESSEL_TABLE_NZ - 2, iZ))
            z_lo = DBLE(iZ) * dz
            alphaZ = (zq - z_lo) / dz
            alphaZ = MAX(0.0D0, MIN(1.0D0, alphaZ))

            f00 = BESSEL_TABLE_PHI(iR+1, iZ+1, ic)
            f10 = BESSEL_TABLE_PHI(iR+2, iZ+1, ic)
            f01 = BESSEL_TABLE_PHI(iR+1, iZ+2, ic)
            f11 = BESSEL_TABLE_PHI(iR+2, iZ+2, ic)

            gr00 = BESSEL_TABLE_DPHI_DR(iR+1, iZ+1, ic)
            gr10 = BESSEL_TABLE_DPHI_DR(iR+2, iZ+1, ic)
            gr01 = BESSEL_TABLE_DPHI_DR(iR+1, iZ+2, ic)
            gr11 = BESSEL_TABLE_DPHI_DR(iR+2, iZ+2, ic)

            gz00 = BESSEL_TABLE_DPHI_DZ(iR+1, iZ+1, ic)
            gz10 = BESSEL_TABLE_DPHI_DZ(iR+2, iZ+1, ic)
            gz01 = BESSEL_TABLE_DPHI_DZ(iR+1, iZ+2, ic)
            gz11 = BESSEL_TABLE_DPHI_DZ(iR+2, iZ+2, ic)

            ! Use bilinear interpolation for robustness. Bicubic-Hermite can
            ! overshoot near sharp gradients and induce non-physical force sign
            ! flips in sparse/edge regions.
            phi_loc = (1.0D0 - alphaR) * (1.0D0 - alphaZ) * f00 + &
                      alphaR * (1.0D0 - alphaZ) * f10 + &
                      (1.0D0 - alphaR) * alphaZ * f01 + &
                      alphaR * alphaZ * f11

            dphi_dR_loc = (1.0D0 - alphaR) * (1.0D0 - alphaZ) * gr00 + &
                          alphaR * (1.0D0 - alphaZ) * gr10 + &
                          (1.0D0 - alphaR) * alphaZ * gr01 + &
                          alphaR * alphaZ * gr11

            dphi_dz_abs = (1.0D0 - alphaR) * (1.0D0 - alphaZ) * gz00 + &
                          alphaR * (1.0D0 - alphaZ) * gz10 + &
                          (1.0D0 - alphaR) * alphaZ * gz01 + &
                          alphaR * alphaZ * gz11

            phi(i) = phi_loc
            aR_val = -dphi_dR_loc
            az(i) = -signz * dphi_dz_abs

            IF (R > eps) THEN
                ax(i) = aR_val * x(i) / R
                ay(i) = aR_val * y(i) / R
            ELSE
                ax(i) = 0.0D0
                ay(i) = 0.0D0
            END IF

            ! Far-field stabilization: smoothly blend to the monopole closure
            ! built from the table-estimated enclosed mass. This guarantees
            ! attractive asymptotics and damps transform truncation oscillations
            ! at large radii.
            scale_ref = MAX(BESSEL_COMPONENT_R_SCALE(ic), BESSEL_COMPONENT_Z_SCALE(ic))
            blend0 = 8.0D0 * scale_ref
            blend1 = 20.0D0 * scale_ref

            r_sph = SQRT(R*R + z(i)*z(i))
            mass_est = MAX(BESSEL_TABLE_META(4, ic), 0.0D0)
            IF (mass_est > 0.0D0 .AND. r_sph > eps .AND. blend1 > blend0) THEN
                t = (r_sph - blend0) / (blend1 - blend0)
                t = MAX(0.0D0, MIN(1.0D0, t))
                w_blend = t*t*(3.0D0 - 2.0D0*t)

                r3 = MAX(r_sph**3, eps)
                phi_mono = -BESSEL_G * mass_est / r_sph

                ax(i) = (1.0D0 - w_blend) * ax(i) + w_blend * (-BESSEL_G * mass_est * x(i) / r3)
                ay(i) = (1.0D0 - w_blend) * ay(i) + w_blend * (-BESSEL_G * mass_est * y(i) / r3)
                az(i) = (1.0D0 - w_blend) * az(i) + w_blend * (-BESSEL_G * mass_est * z(i) / r3)
                phi(i) = (1.0D0 - w_blend) * phi(i) + w_blend * phi_mono
            END IF
        END DO
    END SUBROUTINE component_force_potential

    SUBROUTINE build_table(component_index, density_model, params)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: component_index
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        PROCEDURE(axisymmetric_density_model) :: density_model
        REAL*8, PARAMETER :: pi = 3.14159265358979323846D0
        INTEGER :: iR, iZ, j, m
        REAL*8 :: R, z_val, z1, z2
        REAL*8 :: logR_min, logR_max, dlogR, z_max, dz
        REAL*8 :: sum_phi, sum_dphi_dR, sum_dphi_dz, j0v, j1v
        REAL*8 :: kval, t, mu, kmax, k_scale
        REAL*8 :: r_scale, z_scale, inner1, inner2, d_inner1, d_inner2
        REAL*8 :: e_abs1, e_abs2, e_even1, e_even2, sgn1, sgn2
        REAL*8 :: dr_loc, zint, mass_est, mass_table_est
        REAL*8 :: sigma0, hR_model
        REAL*8, DIMENSION(BESSEL_TABLE_NR) :: r_nodes, int_nodes, zint_nodes, rzint_nodes
        REAL*8, DIMENSION(BESSEL_TABLE_NZ) :: z_nodes, x_line, y_line, rho_line
        REAL*8, DIMENSION(BESSEL_TABLE_NR, BESSEL_TABLE_NZ) :: rho_nodes
        REAL*8, DIMENSION(NK_BUILD, BESSEL_TABLE_NZ) :: rho_kz, zconv_kz, zconv_dz_kz
        REAL*8, DIMENSION(NK_BUILD) :: k_grid, quad_w, mu_q, w_q

        ! Use per-component domain scales to support mixed-scale multi-component models.
        r_scale = MAX(BESSEL_COMPONENT_R_SCALE(component_index), 1.0D-8)
        z_scale = MAX(BESSEL_COMPONENT_Z_SCALE(component_index), 1.0D-8)

        logR_min = LOG(1.0D-3 * r_scale)
        logR_max = LOG(2.0D2 * r_scale)
        dlogR = (logR_max - logR_min) / DBLE(BESSEL_TABLE_NR - 1)
        z_max = 2.0D2 * z_scale
        dz = z_max / DBLE(BESSEL_TABLE_NZ - 1)

        DO iZ = 1, BESSEL_TABLE_NZ
            z_nodes(iZ) = DBLE(iZ - 1) * dz
        END DO

        x_line = 0.0D0
        y_line = 0.0D0

        DO iR = 1, BESSEL_TABLE_NR
            r_nodes(iR) = EXP(logR_min + DBLE(iR - 1) * dlogR)
            x_line = r_nodes(iR)
            CALL density_model(params, BESSEL_TABLE_NZ, x_line, y_line, z_nodes, rho_line)
            rho_nodes(iR, :) = rho_line(:)
        END DO

        BESSEL_TABLE_META(1, component_index) = logR_min
        BESSEL_TABLE_META(2, component_index) = dlogR
        BESSEL_TABLE_META(3, component_index) = dz
        BESSEL_TABLE_META(4, component_index) = 0.0D0

        CALL gauss_legendre_nodes_weights(NK_BUILD, mu_q, w_q)
        ! Bounded k quadrature is more stable in practice than [0,inf) rational mapping
        ! for oscillatory midplane integrals in table construction.
        k_scale = MAX(1.0D0 / r_scale, 1.0D0 / z_scale)
        kmax = 40.0D0 * k_scale
        DO j = 1, NK_BUILD
            mu = mu_q(j)
            t = 0.5D0 * (mu + 1.0D0)
            kval = t * kmax
            k_grid(j) = kval
            quad_w(j) = 0.5D0 * w_q(j) * kmax
        END DO

        ! Radial Hankel transform of rho(R,z) for each z-slice: rho_kz(k,z)
        DO j = 1, NK_BUILD
            kval = k_grid(j)
            DO iZ = 1, BESSEL_TABLE_NZ
                DO iR = 1, BESSEL_TABLE_NR
                    int_nodes(iR) = r_nodes(iR) * bessel_j0_scalar(kval * r_nodes(iR)) * rho_nodes(iR, iZ)
                END DO

                rho_kz(j, iZ) = 0.0D0
                DO iR = 1, BESSEL_TABLE_NR - 1
                    dr_loc = r_nodes(iR+1) - r_nodes(iR)
                    rho_kz(j, iZ) = rho_kz(j, iZ) + 0.5D0 * (int_nodes(iR) + int_nodes(iR+1)) * dr_loc
                END DO
            END DO

            ! Vertical convolution for even-density extension:
            ! int_0^inf rho_k(z') [exp(-k|z-z'|) + exp(-k(z+z'))] dz'
            DO iZ = 1, BESSEL_TABLE_NZ
                z_val = z_nodes(iZ)
                zconv_kz(j, iZ) = 0.0D0
                zconv_dz_kz(j, iZ) = 0.0D0
                DO m = 1, BESSEL_TABLE_NZ - 1
                    z1 = z_nodes(m)
                    z2 = z_nodes(m+1)

                    e_abs1 = EXP(-kval * ABS(z_val - z1))
                    e_abs2 = EXP(-kval * ABS(z_val - z2))
                    e_even1 = EXP(-kval * (z_val + z1))
                    e_even2 = EXP(-kval * (z_val + z2))

                    inner1 = rho_kz(j, m) * (e_abs1 + e_even1)
                    inner2 = rho_kz(j, m+1) * (e_abs2 + e_even2)

                    IF (z_val > z1) THEN
                        sgn1 = 1.0D0
                    ELSE IF (z_val < z1) THEN
                        sgn1 = -1.0D0
                    ELSE
                        sgn1 = 0.0D0
                    END IF
                    IF (z_val > z2) THEN
                        sgn2 = 1.0D0
                    ELSE IF (z_val < z2) THEN
                        sgn2 = -1.0D0
                    ELSE
                        sgn2 = 0.0D0
                    END IF

                    d_inner1 = rho_kz(j, m) * (-kval * sgn1 * e_abs1 - kval * e_even1)
                    d_inner2 = rho_kz(j, m+1) * (-kval * sgn2 * e_abs2 - kval * e_even2)

                    zconv_kz(j, iZ) = zconv_kz(j, iZ) + 0.5D0 * (inner1 + inner2) * (z2 - z1)
                    zconv_dz_kz(j, iZ) = zconv_dz_kz(j, iZ) + 0.5D0 * (d_inner1 + d_inner2) * (z2 - z1)
                END DO
            END DO
        END DO

        DO iR = 1, BESSEL_TABLE_NR
            R = EXP(logR_min + DBLE(iR - 1) * dlogR)
            DO iZ = 1, BESSEL_TABLE_NZ
                sum_phi = 0.0D0
                sum_dphi_dR = 0.0D0
                sum_dphi_dz = 0.0D0
                DO j = 1, NK_BUILD
                    kval = k_grid(j)
                    j0v = bessel_j0_scalar(kval * R)
                    j1v = bessel_j1_scalar(kval * R)

                    sum_phi = sum_phi + quad_w(j) * j0v * zconv_kz(j, iZ)
                    sum_dphi_dR = sum_dphi_dR + quad_w(j) * kval * j1v * zconv_kz(j, iZ)
                    sum_dphi_dz = sum_dphi_dz + quad_w(j) * j0v * zconv_dz_kz(j, iZ)
                END DO

                BESSEL_TABLE_PHI(iR, iZ, component_index) = -2.0D0 * pi * BESSEL_G * sum_phi
                BESSEL_TABLE_DPHI_DR(iR, iZ, component_index) = 2.0D0 * pi * BESSEL_G * sum_dphi_dR
                BESSEL_TABLE_DPHI_DZ(iR, iZ, component_index) = -2.0D0 * pi * BESSEL_G * sum_dphi_dz
            END DO
        END DO

        ! Estimate total mass from sampled density table using even-z symmetry:
        ! M = 4*pi * int_0^inf [ R * int_0^inf rho(R,z) dz ] dR
        ! Use trapezoidal integration consistently in both z and R.
        DO iR = 1, BESSEL_TABLE_NR
            zint = 0.0D0
            DO iZ = 1, BESSEL_TABLE_NZ - 1
                zint = zint + 0.5D0 * (rho_nodes(iR, iZ) + rho_nodes(iR, iZ+1)) * (z_nodes(iZ+1) - z_nodes(iZ))
            END DO
            zint_nodes(iR) = zint
            rzint_nodes(iR) = r_nodes(iR) * zint_nodes(iR)
        END DO

        mass_table_est = 0.0D0
        DO iR = 1, BESSEL_TABLE_NR - 1
            mass_table_est = mass_table_est + 0.5D0 * (rzint_nodes(iR) + rzint_nodes(iR+1)) * (r_nodes(iR+1) - r_nodes(iR))
        END DO
        mass_est = 4.0D0 * pi * mass_table_est

        ! For exponentialdisk params=(Sigma0, hR, hZ), enforce exact total-mass
        ! normalization for far-field monopole closure.
        IF (SIZE(params) >= 3) THEN
            sigma0 = params(1)
            hR_model = params(2)
            IF (sigma0 > 0.0D0 .AND. hR_model > 0.0D0) THEN
                mass_est = 2.0D0 * pi * sigma0 * hR_model * hR_model
            END IF
        END IF

        BESSEL_TABLE_META(4, component_index) = MAX(mass_est, 0.0D0)

        DO iR = 1, BESSEL_TABLE_NR
            DO iZ = 1, BESSEL_TABLE_NZ
                IF (iZ == 1) THEN
                    BESSEL_TABLE_D2PHI_DRDZ(iR, iZ, component_index) = &
                        (BESSEL_TABLE_DPHI_DR(iR, iZ+1, component_index) - &
                         BESSEL_TABLE_DPHI_DR(iR, iZ, component_index)) / dz
                ELSE IF (iZ == BESSEL_TABLE_NZ) THEN
                    BESSEL_TABLE_D2PHI_DRDZ(iR, iZ, component_index) = &
                        (BESSEL_TABLE_DPHI_DR(iR, iZ, component_index) - &
                         BESSEL_TABLE_DPHI_DR(iR, iZ-1, component_index)) / dz
                ELSE
                    BESSEL_TABLE_D2PHI_DRDZ(iR, iZ, component_index) = &
                        (BESSEL_TABLE_DPHI_DR(iR, iZ+1, component_index) - &
                         BESSEL_TABLE_DPHI_DR(iR, iZ-1, component_index)) / (2.0D0 * dz)
                END IF
            END DO
        END DO
    END SUBROUTINE build_table

END MODULE besselbfe
