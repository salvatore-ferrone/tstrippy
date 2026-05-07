MODULE sphericalharmonicsbfe
    USE mathutils, ONLY: linear_interp_scalar, legendre_axisymmetric_basis, legendre_p_all_axisymmetric, gauss_legendre_nodes_weights
    IMPLICIT NONE

    ABSTRACT INTERFACE
        SUBROUTINE axisymmetric_density_model(params, n, x, y, z, rho)
            REAL*8, INTENT(IN), DIMENSION(:) :: params
            INTEGER, INTENT(IN) :: n
            REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
            REAL*8, INTENT(OUT), DIMENSION(n) :: rho
        END SUBROUTINE axisymmetric_density_model
    END INTERFACE

    REAL*8, PRIVATE :: G_SHBFE_DEFAULT = 4.30091727D-6
    REAL*8, PRIVATE :: G_SHBFE = -1.0D0
    LOGICAL, PRIVATE :: SHBFE_G_IS_SET = .FALSE.

    REAL*8, DIMENSION(:),   ALLOCATABLE, PUBLIC :: BASIS_R_GRID
    REAL*8, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: BASIS_RHO_L_GRID
    REAL*8, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: BASIS_PHI_L_GRID
    REAL*8, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: BASIS_DPHI_L_DR_GRID
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: BASIS_PHI_L_COMPONENT_GRID
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: BASIS_DPHIDR_COMP_GRID

    LOGICAL, PUBLIC :: BASIS_GRID_SET = .FALSE.
    LOGICAL, PUBLIC :: BASIS_EXPANSION_INITIALIZED = .FALSE.
    INTEGER, PUBLIC :: BASIS_LMAX = -1
    INTEGER, PUBLIC :: BASIS_NR = -1
    INTEGER, PUBLIC :: BASIS_NCOMP = 0

CONTAINS

    SUBROUTINE setsphericalharmonicbasisgravityconstant(g_in)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: g_in

        IF (g_in <= 0.0D0) THEN
            WRITE(*,'(A)') "WARNING: setsphericalharmonicbasisgravityconstant: G must be positive"
            RETURN
        END IF
        G_SHBFE = g_in
        SHBFE_G_IS_SET = .TRUE.
    END SUBROUTINE setsphericalharmonicbasisgravityconstant

    SUBROUTINE clearsphericalharmonicbasis()
        IMPLICIT NONE

        IF (ALLOCATED(BASIS_R_GRID)) DEALLOCATE(BASIS_R_GRID)
        IF (ALLOCATED(BASIS_RHO_L_GRID)) DEALLOCATE(BASIS_RHO_L_GRID)
        IF (ALLOCATED(BASIS_PHI_L_GRID)) DEALLOCATE(BASIS_PHI_L_GRID)
        IF (ALLOCATED(BASIS_DPHI_L_DR_GRID)) DEALLOCATE(BASIS_DPHI_L_DR_GRID)
        IF (ALLOCATED(BASIS_PHI_L_COMPONENT_GRID)) DEALLOCATE(BASIS_PHI_L_COMPONENT_GRID)
        IF (ALLOCATED(BASIS_DPHIDR_COMP_GRID)) DEALLOCATE(BASIS_DPHIDR_COMP_GRID)

        BASIS_GRID_SET = .FALSE.
        BASIS_EXPANSION_INITIALIZED = .FALSE.
        BASIS_LMAX = -1
        BASIS_NR = -1
        BASIS_NCOMP = 0
        G_SHBFE = -1.0D0
        SHBFE_G_IS_SET = .FALSE.
    END SUBROUTINE clearsphericalharmonicbasis

    SUBROUTINE initsphericalharmonicbasis(lmax, nr, r_grid)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: lmax, nr
        REAL*8, INTENT(IN), DIMENSION(nr) :: r_grid

        CALL clearsphericalharmonicbasis()

        ALLOCATE(BASIS_R_GRID(nr))
        ALLOCATE(BASIS_RHO_L_GRID(0:lmax, nr))
        ALLOCATE(BASIS_PHI_L_GRID(0:lmax, nr))
        ALLOCATE(BASIS_DPHI_L_DR_GRID(0:lmax, nr))

        BASIS_RHO_L_GRID = 0.0D0
        BASIS_PHI_L_GRID = 0.0D0
        BASIS_DPHI_L_DR_GRID = 0.0D0

        BASIS_LMAX = lmax
        BASIS_NR = nr
        BASIS_NCOMP = 0
        BASIS_R_GRID = r_grid
        BASIS_GRID_SET = .TRUE.
        BASIS_EXPANSION_INITIALIZED = .FALSE.

        IF (.NOT. SHBFE_G_IS_SET) THEN
            CALL setsphericalharmonicbasisgravityconstant(G_SHBFE_DEFAULT)
        END IF
    END SUBROUTINE initsphericalharmonicbasis

    SUBROUTINE defaultinitsphericalharmonicbasis()
        IMPLICIT NONE
        INTEGER, PARAMETER :: default_lmax = 20
        INTEGER, PARAMETER :: default_nr = 100
        REAL*8, PARAMETER :: default_rmin = 1.0D-4
        REAL*8, PARAMETER :: default_rmax = 1.0D3
        REAL*8 :: r_grid(default_nr)
        REAL*8 :: log_rmin, log_rmax, dlog_r
        INTEGER :: i

        log_rmin = LOG(default_rmin)
        log_rmax = LOG(default_rmax)
        dlog_r = (log_rmax - log_rmin) / DBLE(default_nr - 1)
        DO i = 1, default_nr
            r_grid(i) = EXP(log_rmin + (i-1) * dlog_r)
        END DO

        CALL initsphericalharmonicbasis(default_lmax, default_nr, r_grid)
    END SUBROUTINE defaultinitsphericalharmonicbasis

    SUBROUTINE initsphericalharmoniccomponentphi(ncomp)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ncomp

        IF (.NOT. BASIS_GRID_SET) THEN
            WRITE(*,'(A)') "WARNING: initsphericalharmoniccomponentphi: call initsphericalharmonicbasis first"
            RETURN
        END IF
        IF (ncomp < 1) THEN
            IF (ALLOCATED(BASIS_PHI_L_COMPONENT_GRID)) DEALLOCATE(BASIS_PHI_L_COMPONENT_GRID)
            IF (ALLOCATED(BASIS_DPHIDR_COMP_GRID)) DEALLOCATE(BASIS_DPHIDR_COMP_GRID)
            BASIS_NCOMP = 0
            RETURN
        END IF

        IF (ALLOCATED(BASIS_PHI_L_COMPONENT_GRID)) DEALLOCATE(BASIS_PHI_L_COMPONENT_GRID)
        IF (ALLOCATED(BASIS_DPHIDR_COMP_GRID)) DEALLOCATE(BASIS_DPHIDR_COMP_GRID)
        ALLOCATE(BASIS_PHI_L_COMPONENT_GRID(0:BASIS_LMAX, BASIS_NR, ncomp))
        ALLOCATE(BASIS_DPHIDR_COMP_GRID(0:BASIS_LMAX, BASIS_NR, ncomp))
        BASIS_PHI_L_COMPONENT_GRID = 0.0D0
        BASIS_DPHIDR_COMP_GRID = 0.0D0
        BASIS_NCOMP = ncomp
    END SUBROUTINE initsphericalharmoniccomponentphi

    SUBROUTINE storesphericalharmoniccomponentphi(icomp)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: icomp

        IF (.NOT. BASIS_GRID_SET) THEN
            WRITE(*,'(A)') "WARNING: storesphericalharmoniccomponentphi: call initsphericalharmonicbasis first"
            RETURN
        END IF
        IF (.NOT. ALLOCATED(BASIS_PHI_L_COMPONENT_GRID)) THEN
            WRITE(*,'(A)') "WARNING: storesphericalharmoniccomponentphi: call initsphericalharmoniccomponentphi first"
            RETURN
        END IF
        IF (icomp < 1 .OR. icomp > BASIS_NCOMP) THEN
            WRITE(*,'(A)') "WARNING: storesphericalharmoniccomponentphi: component index out of range"
            RETURN
        END IF

        BASIS_PHI_L_COMPONENT_GRID(:,:,icomp) = BASIS_PHI_L_GRID(:,:)
        BASIS_DPHIDR_COMP_GRID(:,:,icomp) = BASIS_DPHI_L_DR_GRID(:,:)
    END SUBROUTINE storesphericalharmoniccomponentphi

    SUBROUTINE loadsphericalharmoniccomponentphi(icomp)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: icomp

        IF (.NOT. BASIS_GRID_SET) THEN
            WRITE(*,'(A)') "WARNING: loadsphericalharmoniccomponentphi: call initsphericalharmonicbasis first"
            RETURN
        END IF
        IF (.NOT. ALLOCATED(BASIS_PHI_L_COMPONENT_GRID)) THEN
            WRITE(*,'(A)') "WARNING: loadsphericalharmoniccomponentphi: call initsphericalharmoniccomponentphi first"
            RETURN
        END IF
        IF (icomp < 1 .OR. icomp > BASIS_NCOMP) THEN
            WRITE(*,'(A)') "WARNING: loadsphericalharmoniccomponentphi: component index out of range"
            RETURN
        END IF

        BASIS_PHI_L_GRID(:,:) = BASIS_PHI_L_COMPONENT_GRID(:,:,icomp)
        BASIS_DPHI_L_DR_GRID(:,:) = BASIS_DPHIDR_COMP_GRID(:,:,icomp)
    END SUBROUTINE loadsphericalharmoniccomponentphi

    SUBROUTINE project_axisym_density_generic(params, density_model)
        ! Generic Legendre projection for any axisymmetric density model.
        IMPLICIT NONE
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        PROCEDURE(axisymmetric_density_model) :: density_model
        INTEGER :: n_mu, i_r, l, k
        REAL*8, ALLOCATABLE :: mu_q(:), w_q(:), p(:)
        REAL*8 :: mu, factor
        REAL*8 :: x1(1), y1(1), z1(1), rho1(1), r, rcyl

        IF (.NOT. BASIS_GRID_SET) THEN
            WRITE(*,'(A)') "WARNING: project_axisym_density_generic: call initsphericalharmonicbasis first"
            RETURN
        END IF

        n_mu = MAX(4 * (BASIS_LMAX + 1), 40)
        ALLOCATE(mu_q(n_mu), w_q(n_mu), p(0:BASIS_LMAX))
        CALL gauss_legendre_nodes_weights(n_mu, mu_q, w_q)

        BASIS_RHO_L_GRID = 0.0D0
        DO i_r = 1, BASIS_NR
            r = BASIS_R_GRID(i_r)
            DO k = 1, n_mu
                mu = mu_q(k)
                rcyl = r * SQRT(MAX(1.0D0 - mu*mu, 0.0D0))
                x1(1) = rcyl
                y1(1) = 0.0D0
                z1(1) = r * mu
                CALL density_model(params, 1, x1, y1, z1, rho1)

                CALL legendre_p_all_axisymmetric(BASIS_LMAX, mu, p)
                DO l = 0, BASIS_LMAX, 2
                    factor = (2*l + 1) * 0.5D0 * w_q(k)
                    BASIS_RHO_L_GRID(l, i_r) = BASIS_RHO_L_GRID(l, i_r) + factor * rho1(1) * p(l)
                END DO
            END DO
        END DO

        DEALLOCATE(mu_q, w_q, p)
    END SUBROUTINE project_axisym_density_generic

    SUBROUTINE compute_phi_tables_from_rho()
        IMPLICIT NONE
        REAL*8, PARAMETER :: pi_phi = 3.14159265358979323846D0
        INTEGER :: l, i
        REAL*8 :: prefactor, r_lo, r_hi, f_lo, f_hi, dr
        REAL*8, ALLOCATABLE :: i_less(:), i_greater(:)

        IF (.NOT. BASIS_GRID_SET) THEN
            WRITE(*,'(A)') "WARNING: compute_phi_tables_from_rho: basis not initialized"
            RETURN
        END IF

        ALLOCATE(i_less(BASIS_NR), i_greater(BASIS_NR))

        DO l = 0, BASIS_LMAX, 2
            prefactor = -4.0D0 * pi_phi * G_SHBFE / DBLE(2*l + 1)

            i_less(1) = 0.0D0
            DO i = 2, BASIS_NR
                r_lo = BASIS_R_GRID(i-1)
                r_hi = BASIS_R_GRID(i)
                f_lo = r_lo**(l+2) * BASIS_RHO_L_GRID(l, i-1)
                f_hi = r_hi**(l+2) * BASIS_RHO_L_GRID(l, i)
                i_less(i) = i_less(i-1) + 0.5D0 * (r_hi - r_lo) * (f_lo + f_hi)
            END DO

            i_greater(BASIS_NR) = 0.0D0
            DO i = BASIS_NR-1, 1, -1
                r_lo = BASIS_R_GRID(i)
                r_hi = BASIS_R_GRID(i+1)
                f_lo = r_lo**(1-l) * BASIS_RHO_L_GRID(l, i)
                f_hi = r_hi**(1-l) * BASIS_RHO_L_GRID(l, i+1)
                dr = r_hi - r_lo
                i_greater(i) = i_greater(i+1) + 0.5D0 * dr * (f_lo + f_hi)
            END DO

            DO i = 1, BASIS_NR
                BASIS_PHI_L_GRID(l, i) = prefactor * ( &
                    BASIS_R_GRID(i)**(-(l+1)) * i_less(i) + &
                    BASIS_R_GRID(i)**l * i_greater(i) )
                BASIS_DPHI_L_DR_GRID(l, i) = prefactor * ( &
                    -(l+1) * BASIS_R_GRID(i)**(-(l+2)) * i_less(i) + &
                    DBLE(l) * BASIS_R_GRID(i)**(l-1) * i_greater(i) )
            END DO
        END DO

        DEALLOCATE(i_less, i_greater)
        BASIS_EXPANSION_INITIALIZED = .TRUE.
    END SUBROUTINE compute_phi_tables_from_rho

    SUBROUTINE sphericalharmonicbasisforce(n, x, y, z, ax, ay, az)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: ax, ay, az
        INTEGER :: i, l, jlo, jhi, jmid
        REAL*8 :: r, r_safe, mu, rcyl2, alpha_r
        REAL*8 :: phi_l_r, dphi_l_dr_r, dphi_dr_v, dphi_dmu_v
        REAL*8 :: p(0:BASIS_LMAX), dp_dmu(0:BASIS_LMAX)
        REAL*8, PARAMETER :: eps_r = 1.0D-30

        DO i = 1, n
            r = SQRT(x(i)**2 + y(i)**2 + z(i)**2)
            r_safe = MAX(r, eps_r)
            mu = z(i) / r_safe
            rcyl2 = x(i)**2 + y(i)**2

            IF (r_safe <= BASIS_R_GRID(1)) THEN
                jlo = 1
                alpha_r = 0.0D0
            ELSE IF (r_safe >= BASIS_R_GRID(BASIS_NR)) THEN
                jlo = BASIS_NR - 1
                alpha_r = 1.0D0
            ELSE
                jlo = 1
                jhi = BASIS_NR
                DO WHILE (jhi - jlo > 1)
                    jmid = (jlo + jhi) / 2
                    IF (BASIS_R_GRID(jmid) <= r_safe) THEN
                        jlo = jmid
                    ELSE
                        jhi = jmid
                    END IF
                END DO
                alpha_r = LOG(r_safe / BASIS_R_GRID(jlo)) / LOG(BASIS_R_GRID(jlo+1) / BASIS_R_GRID(jlo))
            END IF

            CALL legendre_axisymmetric_basis(BASIS_LMAX, mu, p, dp_dmu)

            dphi_dr_v = 0.0D0
            dphi_dmu_v = 0.0D0
            DO l = 0, BASIS_LMAX, 2
                phi_l_r = linear_interp_scalar(BASIS_PHI_L_GRID(l,jlo), BASIS_PHI_L_GRID(l,jlo+1), alpha_r)
                dphi_l_dr_r = linear_interp_scalar(BASIS_DPHI_L_DR_GRID(l,jlo), BASIS_DPHI_L_DR_GRID(l,jlo+1), alpha_r)
                dphi_dr_v = dphi_dr_v + dphi_l_dr_r * p(l)
                dphi_dmu_v = dphi_dmu_v + phi_l_r * dp_dmu(l)
            END DO

            ax(i) = -dphi_dr_v * x(i)/r_safe + dphi_dmu_v * z(i)*x(i)/r_safe**3
            ay(i) = -dphi_dr_v * y(i)/r_safe + dphi_dmu_v * z(i)*y(i)/r_safe**3
            az(i) = -dphi_dr_v * z(i)/r_safe - dphi_dmu_v * rcyl2/r_safe**3
        END DO
    END SUBROUTINE sphericalharmonicbasisforce

    SUBROUTINE sphericalharmonicbasispotential(n, x, y, z, phi_out)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: phi_out
        INTEGER :: i, l, jlo, jhi, jmid
        REAL*8 :: r, r_safe, mu, alpha_r, phi_v, phi_l_r
        REAL*8 :: p(0:BASIS_LMAX), dp_dummy(0:BASIS_LMAX)
        REAL*8, PARAMETER :: eps_r = 1.0D-30

        DO i = 1, n
            r = SQRT(x(i)**2 + y(i)**2 + z(i)**2)
            r_safe = MAX(r, eps_r)
            mu = z(i) / r_safe

            IF (r_safe <= BASIS_R_GRID(1)) THEN
                jlo = 1
                alpha_r = 0.0D0
            ELSE IF (r_safe >= BASIS_R_GRID(BASIS_NR)) THEN
                jlo = BASIS_NR - 1
                alpha_r = 1.0D0
            ELSE
                jlo = 1
                jhi = BASIS_NR
                DO WHILE (jhi - jlo > 1)
                    jmid = (jlo + jhi) / 2
                    IF (BASIS_R_GRID(jmid) <= r_safe) THEN
                        jlo = jmid
                    ELSE
                        jhi = jmid
                    END IF
                END DO
                alpha_r = LOG(r_safe / BASIS_R_GRID(jlo)) / LOG(BASIS_R_GRID(jlo+1) / BASIS_R_GRID(jlo))
            END IF

            CALL legendre_axisymmetric_basis(BASIS_LMAX, mu, p, dp_dummy)

            phi_v = 0.0D0
            DO l = 0, BASIS_LMAX, 2
                phi_l_r = linear_interp_scalar(BASIS_PHI_L_GRID(l,jlo), BASIS_PHI_L_GRID(l,jlo+1), alpha_r)
                phi_v = phi_v + phi_l_r * p(l)
            END DO
            phi_out(i) = phi_v
        END DO
    END SUBROUTINE sphericalharmonicbasispotential

END MODULE sphericalharmonicsbfe