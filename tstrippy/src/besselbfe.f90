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
    INTEGER, PARAMETER, PRIVATE :: BESSEL_META_FIELDS = 3

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

    SUBROUTINE initialize(nr, nz, nk_build_in, r_scale, z_scale)
        ! Explicit initialization with all tunable parameters
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: nr, nz, nk_build_in
        REAL*8, INTENT(IN) :: r_scale, z_scale

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
        IF (r_scale <= 0.0D0) THEN
            WRITE(*,'(A)') "WARNING: initialize: r_scale must be positive"
            RETURN
        END IF
        IF (z_scale <= 0.0D0) THEN
            WRITE(*,'(A)') "WARNING: initialize: z_scale must be positive"
            RETURN
        END IF

        CALL clear()
        BESSEL_TABLE_NR = nr
        BESSEL_TABLE_NZ = nz
        NK_BUILD = nk_build_in
        BESSEL_R_SCALE = r_scale
        BESSEL_Z_SCALE = z_scale
        BESSEL_INITIALIZED = .TRUE.
    END SUBROUTINE initialize

    SUBROUTINE default_initialize()
        ! Default initialization using hardcoded parameter defaults
        IMPLICIT NONE
        CALL initialize(BESSEL_TABLE_NR_DEFAULT, BESSEL_TABLE_NZ_DEFAULT, &
                         NK_BUILD_DEFAULT, BESSEL_R_SCALE_DEFAULT, BESSEL_Z_SCALE_DEFAULT)
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
        ALLOCATE(BESSEL_COMPONENT_READY(BESSEL_NCOMP))

        BESSEL_TABLE_PHI = 0.0D0
        BESSEL_TABLE_DPHI_DR = 0.0D0
        BESSEL_TABLE_DPHI_DZ = 0.0D0
        BESSEL_TABLE_D2PHI_DRDZ = 0.0D0
        BESSEL_TABLE_META = 0.0D0
        BESSEL_COMPONENT_READY = .FALSE.
    END SUBROUTINE allocate_component_tables

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
        REAL*8 :: logR, logR_min, dlogR, dz, logR_max
        REAL*8 :: alphaR, alphaZ, zq
        REAL*8 :: R_lo, R_hi, dR, z_lo
        REAL*8 :: phi_loc, dphi_dR_loc, dphi_dz_abs
        REAL*8 :: f00, f10, f01, f11
        REAL*8 :: fx00, fx10, fx01, fx11
        REAL*8 :: fy00, fy10, fy01, fy11
        REAL*8 :: fxy00, fxy10, fxy01, fxy11

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

            fx00 = BESSEL_TABLE_DPHI_DR(iR+1, iZ+1, ic)
            fx10 = BESSEL_TABLE_DPHI_DR(iR+2, iZ+1, ic)
            fx01 = BESSEL_TABLE_DPHI_DR(iR+1, iZ+2, ic)
            fx11 = BESSEL_TABLE_DPHI_DR(iR+2, iZ+2, ic)

            fy00 = BESSEL_TABLE_DPHI_DZ(iR+1, iZ+1, ic)
            fy10 = BESSEL_TABLE_DPHI_DZ(iR+2, iZ+1, ic)
            fy01 = BESSEL_TABLE_DPHI_DZ(iR+1, iZ+2, ic)
            fy11 = BESSEL_TABLE_DPHI_DZ(iR+2, iZ+2, ic)

            fxy00 = BESSEL_TABLE_D2PHI_DRDZ(iR+1, iZ+1, ic)
            fxy10 = BESSEL_TABLE_D2PHI_DRDZ(iR+2, iZ+1, ic)
            fxy01 = BESSEL_TABLE_D2PHI_DRDZ(iR+1, iZ+2, ic)
            fxy11 = BESSEL_TABLE_D2PHI_DRDZ(iR+2, iZ+2, ic)

            CALL bicubic_hermite_eval_2d(f00, f10, f01, f11, &
                                         fx00, fx10, fx01, fx11, &
                                         fy00, fy10, fy01, fy11, &
                                         fxy00, fxy10, fxy01, fxy11, &
                                         dR, dz, alphaR, alphaZ, &
                                         phi_loc, dphi_dR_loc, dphi_dz_abs)

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
        REAL*8 :: sum_phi, j0v
        REAL*8 :: kval, t, mu, kmax, k_scale
        REAL*8 :: r_scale, z_scale, inner1, inner2
        REAL*8 :: dr_loc
        REAL*8, DIMENSION(BESSEL_TABLE_NR) :: r_nodes, int_nodes
        REAL*8, DIMENSION(BESSEL_TABLE_NZ) :: z_nodes, x_line, y_line, rho_line
        REAL*8, DIMENSION(BESSEL_TABLE_NR, BESSEL_TABLE_NZ) :: rho_nodes
        REAL*8, DIMENSION(NK_BUILD, BESSEL_TABLE_NZ) :: rho_kz, zconv_kz
        REAL*8, DIMENSION(NK_BUILD) :: k_grid, quad_w, mu_q, w_q

        ! Use public module-level domain scales (profile-agnostic)
        r_scale = MAX(BESSEL_R_SCALE, 1.0D-8)
        z_scale = MAX(BESSEL_Z_SCALE, 1.0D-8)

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
                DO m = 1, BESSEL_TABLE_NZ - 1
                    z1 = z_nodes(m)
                    z2 = z_nodes(m+1)

                    inner1 = rho_kz(j, m) * (EXP(-kval * ABS(z_val - z1)) + EXP(-kval * (z_val + z1)))
                    inner2 = rho_kz(j, m+1) * (EXP(-kval * ABS(z_val - z2)) + EXP(-kval * (z_val + z2)))

                    zconv_kz(j, iZ) = zconv_kz(j, iZ) + 0.5D0 * (inner1 + inner2) * (z2 - z1)
                END DO
            END DO
        END DO

        DO iR = 1, BESSEL_TABLE_NR
            R = EXP(logR_min + DBLE(iR - 1) * dlogR)
            DO iZ = 1, BESSEL_TABLE_NZ
                sum_phi = 0.0D0
                DO j = 1, NK_BUILD
                    kval = k_grid(j)
                    j0v = bessel_j0_scalar(kval * R)

                    sum_phi = sum_phi + quad_w(j) * j0v * zconv_kz(j, iZ)
                END DO

                BESSEL_TABLE_PHI(iR, iZ, component_index) = -2.0D0 * pi * BESSEL_G * sum_phi
            END DO
        END DO

        ! Compute dphi/dR from the potential table on the nonuniform R grid.
        DO iZ = 1, BESSEL_TABLE_NZ
            DO iR = 1, BESSEL_TABLE_NR
                IF (iR == 1) THEN
                    BESSEL_TABLE_DPHI_DR(iR, iZ, component_index) = &
                        (BESSEL_TABLE_PHI(iR+1, iZ, component_index) - &
                         BESSEL_TABLE_PHI(iR, iZ, component_index)) / &
                        MAX(r_nodes(iR+1) - r_nodes(iR), 1.0D-16)
                ELSE IF (iR == BESSEL_TABLE_NR) THEN
                    BESSEL_TABLE_DPHI_DR(iR, iZ, component_index) = &
                        (BESSEL_TABLE_PHI(iR, iZ, component_index) - &
                         BESSEL_TABLE_PHI(iR-1, iZ, component_index)) / &
                        MAX(r_nodes(iR) - r_nodes(iR-1), 1.0D-16)
                ELSE
                    BESSEL_TABLE_DPHI_DR(iR, iZ, component_index) = &
                        (BESSEL_TABLE_PHI(iR+1, iZ, component_index) - &
                         BESSEL_TABLE_PHI(iR-1, iZ, component_index)) / &
                        MAX(r_nodes(iR+1) - r_nodes(iR-1), 1.0D-16)
                END IF
            END DO
        END DO

        ! Compute dphi/dz from the potential table for consistency with the
        ! chosen vertical Green-function discretization.
        DO iR = 1, BESSEL_TABLE_NR
            DO iZ = 1, BESSEL_TABLE_NZ
                IF (iZ == 1) THEN
                    BESSEL_TABLE_DPHI_DZ(iR, iZ, component_index) = &
                        (BESSEL_TABLE_PHI(iR, iZ+1, component_index) - &
                         BESSEL_TABLE_PHI(iR, iZ, component_index)) / dz
                ELSE IF (iZ == BESSEL_TABLE_NZ) THEN
                    BESSEL_TABLE_DPHI_DZ(iR, iZ, component_index) = &
                        (BESSEL_TABLE_PHI(iR, iZ, component_index) - &
                         BESSEL_TABLE_PHI(iR, iZ-1, component_index)) / dz
                ELSE
                    BESSEL_TABLE_DPHI_DZ(iR, iZ, component_index) = &
                        (BESSEL_TABLE_PHI(iR, iZ+1, component_index) - &
                         BESSEL_TABLE_PHI(iR, iZ-1, component_index)) / (2.0D0 * dz)
                END IF
            END DO
        END DO

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
