MODULE gravity
    USE mathutils, ONLY: linear_interp_scalar, bicubic_hermite_eval_2d, legendre_axisymmetric_basis, &
                         legendre_p_all_axisymmetric, gauss_legendre_nodes_weights, &
                         bessel_j0_scalar, bessel_j1_scalar
    IMPLICIT NONE

    ! =========================================================================
    ! ORGANIZATION
    ! =========================================================================
    ! This module provides gravitational potential and force evaluators
    ! organized into logical sections:
    !
    ! 1. STATE SUBROUTINES — Lifecycle API (cleargravity, addgravitycomponent, etc.)
    ! 2. ANALYTICAL POTENTIALS — Closed-form models (Plummer, Hernquist, etc.)
    ! 3. DENSITY-ONLY POTENTIALS — Legendre expansion models (exponential halo, etc.)
    ! 4. BFE INFRASTRUCTURE — Composite basis expansion setup and evaluation
    ! 5. SPHERICAL HARMONICS — Legendre projection and evaluation primitives
    ! 6. BESSEL/TABLE INFRASTRUCTURE — Disk table construction and evaluation
    !
    ! All public symbols are listed in module state or at section starts.
    ! Private symbols are internal implementation details and marked PRIVATE.
    ! =========================================================================

    ABSTRACT INTERFACE
        SUBROUTINE bessel_component_evaluator_interface(params, N, x, y, z, ax, ay, az, phi)
            REAL*8, INTENT(IN), DIMENSION(*) :: params
            INTEGER, INTENT(IN) :: N
            REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
            REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi
        END SUBROUTINE bessel_component_evaluator_interface
    END INTERFACE

    ! Axisymmetric basis-expansion tables (l-mode by radial grid index).
    ! BASIS_GRID_SET: r_grid allocated, arrays ready for projection.
    ! BASIS_EXPANSION_INITIALIZED: projection + potential tables computed, ready to evaluate.
    LOGICAL, PUBLIC  :: BASIS_GRID_SET              = .FALSE.
    LOGICAL, PUBLIC  :: BASIS_EXPANSION_INITIALIZED = .FALSE.
    INTEGER, PUBLIC  :: BASIS_LMAX = -1
    INTEGER, PUBLIC  :: BASIS_NR   = -1
    REAL*8,  PUBLIC  :: BASIS_G    = -1.0D0

    REAL*8, DIMENSION(:),   ALLOCATABLE, PUBLIC :: BASIS_R_GRID
    REAL*8, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: BASIS_RHO_L_GRID
    REAL*8, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: BASIS_PHI_L_GRID
    REAL*8, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: BASIS_DPHI_L_DR_GRID

    ! Composite BFE manager state.
    INTEGER, PARAMETER, PUBLIC :: BFE_KIND_NONE = 0
    INTEGER, PARAMETER, PUBLIC :: BFE_KIND_LEGENDRE = 1
    INTEGER, PARAMETER, PUBLIC :: BFE_KIND_BESSEL_DISK = 2    ! direct quadrature (reference only)
    INTEGER, PARAMETER, PUBLIC :: BFE_KIND_DISK_TABLE = 3     ! precomputed cylindrical table (production)
    INTEGER, PARAMETER, PUBLIC :: COMPOSITE_BESSEL_MAX_PARAMS = 16

    ! Disk table resolution (shared across all disk components).
    INTEGER, PARAMETER, PUBLIC :: COMPOSITE_DISK_TABLE_NR = 256  ! log-spaced R nodes
    INTEGER, PARAMETER, PUBLIC :: COMPOSITE_DISK_TABLE_NZ = 128  ! linear z nodes (z >= 0)

    LOGICAL, PUBLIC :: COMPOSITE_BASIS_GRID_SET = .FALSE.
    LOGICAL, PUBLIC :: COMPOSITE_BASIS_FINALIZED = .FALSE.
    INTEGER, PUBLIC :: COMPOSITE_NCOMP = 0
    INTEGER, PUBLIC :: COMPOSITE_LMAX = -1
    INTEGER, PUBLIC :: COMPOSITE_NR = -1

    INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: COMPOSITE_KIND
    LOGICAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: COMPOSITE_READY
    INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: COMPOSITE_BESSEL_NPARAMS

    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: COMPOSITE_R_GRID
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: COMPOSITE_RHO_L_GRID
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: COMPOSITE_PHI_L_GRID
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: COMPOSITE_DPHI_L_DR_GRID

    REAL*8, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: COMPOSITE_BESSEL_PARAMS

    ! Per-component disk-table metadata: row 1 = logR_min, row 2 = dlogR, row 3 = dz.
    REAL*8, DIMENSION(:,:),   ALLOCATABLE, PUBLIC :: COMPOSITE_DISK_TABLE_META   ! (3, ncomp)
    ! Precomputed cylindrical (R, z) tables per disk component.
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: COMPOSITE_DISK_TABLE_PHI        ! (nR, nZ, ncomp)
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: COMPOSITE_DISK_TABLE_DPHI_DR    ! (nR, nZ, ncomp)
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: COMPOSITE_DISK_TABLE_DPHI_DZ    ! (nR, nZ, ncomp), z>=0
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: COMPOSITE_DISK_TABLE_D2PHI_DRDZ ! (nR, nZ, ncomp), z>=0

    ! -----------------------------------------------------------------------
    ! Gravity lifecycle module state (unified API)
    ! -----------------------------------------------------------------------
    REAL*8,  PARAMETER, PUBLIC :: GRAVITY_G_DEFAULT = 4.30091727D-6  ! kpc (km/s)^2 / Msun
    INTEGER, PARAMETER, PUBLIC :: GRAVITY_MAX_NCOMP  = 16
    INTEGER, PARAMETER, PUBLIC :: GRAVITY_MAX_PARAMS = 16

    INTEGER, PARAMETER, PUBLIC :: GRAVITY_KIND_NONE         = 0
    INTEGER, PARAMETER, PUBLIC :: GRAVITY_KIND_PLUMMER       = 10
    INTEGER, PARAMETER, PUBLIC :: GRAVITY_KIND_HERNQUIST     = 11
    INTEGER, PARAMETER, PUBLIC :: GRAVITY_KIND_MIYAMOTONAGAI = 12
    INTEGER, PARAMETER, PUBLIC :: GRAVITY_KIND_LONGMURALIBAR = 13

    REAL*8,  PUBLIC :: GRAVITY_G            = GRAVITY_G_DEFAULT
    LOGICAL, PUBLIC :: GRAVITY_G_IS_DEFAULT = .TRUE.
    LOGICAL, PUBLIC :: GRAVITY_FINALIZED    = .FALSE.
    INTEGER, PUBLIC :: GRAVITY_NCOMP        = 0
    INTEGER, DIMENSION(GRAVITY_MAX_NCOMP),                    PUBLIC :: GRAVITY_KIND   = 0
    REAL*8,  DIMENSION(GRAVITY_MAX_PARAMS, GRAVITY_MAX_NCOMP), PUBLIC :: GRAVITY_PARAMS = 0.0D0

    ! Private internal routines (not exposed to Python/caller)
    ! NOTE: No backward compatibility aliases. Use canonical names only.
    PRIVATE :: default_init_basis_expansion
    PRIVATE :: project_exponential_oblate_halo
    PRIVATE :: project_ibata2024halo
    PRIVATE :: compute_phi_tables_from_rho
    PRIVATE :: compute_phi_tables_from_rho_component
    PRIVATE :: axisymmetricbasisexpansion_eval
    PRIVATE :: axisymmetricbasisexpansion_eval_component
    PRIVATE :: bessel_eval_component
    PRIVATE :: exponential_disk_bessel_eval_component
    PRIVATE :: build_exponential_disk_table
    PRIVATE :: disk_table_eval_component

    CONTAINS

    ! =======================================================================
    ! STATE SUBROUTINES (Lifecycle API)
    ! cleargravity, setgravityconstant, addgravitycomponent, finalizegravity,
    ! evaluategravityforces, evaluategravitypotential, printgravitystate
    ! =======================================================================

    SUBROUTINE cleargravity()
        IMPLICIT NONE
        GRAVITY_G            = GRAVITY_G_DEFAULT
        GRAVITY_G_IS_DEFAULT = .TRUE.
        GRAVITY_FINALIZED    = .FALSE.
        GRAVITY_NCOMP        = 0
        GRAVITY_KIND         = GRAVITY_KIND_NONE
        GRAVITY_PARAMS       = 0.0D0
    END SUBROUTINE cleargravity

    SUBROUTINE setgravityconstant(G)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: G
        IF (GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: setgravityconstant: cannot change G after finalizegravity"
            RETURN
        END IF
        IF (G <= 0.0D0) THEN
            WRITE(*,'(A)') "WARNING: setgravityconstant: G must be positive"
            RETURN
        END IF
        GRAVITY_G            = G
        GRAVITY_G_IS_DEFAULT = .FALSE.
    END SUBROUTINE setgravityconstant

    SUBROUTINE addgravitycomponent(model_name, params, nparams)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: model_name
        INTEGER, INTENT(IN) :: nparams
        REAL*8, INTENT(IN), DIMENSION(nparams) :: params
        INTEGER :: kind_code, required_params

        IF (GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: addgravitycomponent: cannot add components after finalizegravity"
            RETURN
        END IF
        IF (GRAVITY_NCOMP >= GRAVITY_MAX_NCOMP) THEN
            WRITE(*,'(A)') "WARNING: addgravitycomponent: maximum number of components reached"
            RETURN
        END IF

        ! Map model name to internal kind code
        SELECT CASE (TRIM(model_name))
        CASE ("plummer")
            kind_code = GRAVITY_KIND_PLUMMER
            required_params = 2  ! M, a
        CASE ("hernquist")
            kind_code = GRAVITY_KIND_HERNQUIST
            required_params = 2  ! M, a
        CASE ("miyamotonagai")
            kind_code = GRAVITY_KIND_MIYAMOTONAGAI
            required_params = 3  ! M, a, b
        CASE ("longmuralibar")
            kind_code = GRAVITY_KIND_LONGMURALIBAR
            required_params = 4  ! M, abar, bbar, cbar
        CASE DEFAULT
            WRITE(*,'(A)') "WARNING: addgravitycomponent: unknown model_name (must be lowercase)"
            RETURN
        END SELECT

        IF (nparams /= required_params) THEN
            WRITE(*,'(A,I0,A,I0)') "WARNING: addgravitycomponent: expected", required_params, &
                " parameters but got", nparams
            RETURN
        END IF

        GRAVITY_NCOMP = GRAVITY_NCOMP + 1
        GRAVITY_KIND(GRAVITY_NCOMP) = kind_code
        GRAVITY_PARAMS(1:nparams, GRAVITY_NCOMP) = params(1:nparams)
    END SUBROUTINE addgravitycomponent

    SUBROUTINE finalizegravity()
        IMPLICIT NONE
        INTEGER :: i

        IF (GRAVITY_NCOMP < 1) THEN
            WRITE(*,'(A)') "WARNING: finalizegravity: no components registered"
            RETURN
        END IF
        DO i = 1, GRAVITY_NCOMP
            IF (GRAVITY_KIND(i) == GRAVITY_KIND_NONE) THEN
                WRITE(*,'(A)') "WARNING: finalizegravity: component slot is uninitialized"
                RETURN
            END IF
        END DO

        IF (GRAVITY_G_IS_DEFAULT) THEN
            WRITE(*,'(A,ES14.8,A)') "finalizegravity: Default G = ", GRAVITY_G, &
                " kpc (km/s)^2 / Msun"
        ELSE
            WRITE(*,'(A,ES14.8,A)') "finalizegravity: G = ", GRAVITY_G, &
                " User inputted value"
        END IF

        GRAVITY_FINALIZED = .TRUE.
    END SUBROUTINE finalizegravity

    SUBROUTINE evaluategravityforces(N, x, y, z, ax, ay, az)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN),  DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az
        REAL*8, DIMENSION(N) :: ax_c, ay_c, az_c
        REAL*8, DIMENSION(N,3) :: force_c
        REAL*8, DIMENSION(2) :: p2
        REAL*8, DIMENSION(3) :: p3
        REAL*8, DIMENSION(4) :: p4
        INTEGER :: i

        IF (.NOT. GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: evaluategravityforces: not finalized yet"
            ax = 0.0D0;  ay = 0.0D0;  az = 0.0D0
            RETURN
        END IF

        ax = 0.0D0;  ay = 0.0D0;  az = 0.0D0

        DO i = 1, GRAVITY_NCOMP
            force_c = 0.0D0
            SELECT CASE (GRAVITY_KIND(i))
            CASE (GRAVITY_KIND_PLUMMER)
                p2 = [GRAVITY_PARAMS(1,i), GRAVITY_PARAMS(2,i)]
                CALL plummerforce(p2, N, x, y, z, force_c)
            CASE (GRAVITY_KIND_HERNQUIST)
                p2 = [GRAVITY_PARAMS(1,i), GRAVITY_PARAMS(2,i)]
                CALL hernquistforce(p2, N, x, y, z, force_c)
            CASE (GRAVITY_KIND_MIYAMOTONAGAI)
                p3 = [GRAVITY_PARAMS(1,i), GRAVITY_PARAMS(2,i), GRAVITY_PARAMS(3,i)]
                CALL miyamotonagaiforce(p3, N, x, y, z, force_c)
            CASE (GRAVITY_KIND_LONGMURALIBAR)
                p4 = [GRAVITY_PARAMS(1,i), GRAVITY_PARAMS(2,i), &
                      GRAVITY_PARAMS(3,i), GRAVITY_PARAMS(4,i)]
                CALL longmuralibarforce(p4, N, x, y, z, force_c)
            CASE DEFAULT
                WRITE(*,'(A)') "WARNING: evaluategravityforces: unknown component kind"
            END SELECT
            ax_c = force_c(:,1)
            ay_c = force_c(:,2)
            az_c = force_c(:,3)
            ax = ax + ax_c
            ay = ay + ay_c
            az = az + az_c
        END DO
    END SUBROUTINE evaluategravityforces

    SUBROUTINE evaluategravitypotential(N, x, y, z, phi)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN),  DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: phi
        REAL*8, DIMENSION(N) :: phi_c
        REAL*8, DIMENSION(2) :: p2
        REAL*8, DIMENSION(3) :: p3
        REAL*8, DIMENSION(4) :: p4
        INTEGER :: i

        IF (.NOT. GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: evaluategravitypotential: not finalized yet"
            phi = 0.0D0
            RETURN
        END IF

        phi = 0.0D0

        DO i = 1, GRAVITY_NCOMP
            phi_c = 0.0D0
            SELECT CASE (GRAVITY_KIND(i))
            CASE (GRAVITY_KIND_PLUMMER)
                p2 = [GRAVITY_PARAMS(1,i), GRAVITY_PARAMS(2,i)]
                CALL plummerpotential(p2, N, x, y, z, phi_c)
            CASE (GRAVITY_KIND_HERNQUIST)
                p2 = [GRAVITY_PARAMS(1,i), GRAVITY_PARAMS(2,i)]
                CALL hernquistpotential(p2, N, x, y, z, phi_c)
            CASE (GRAVITY_KIND_MIYAMOTONAGAI)
                p3 = [GRAVITY_PARAMS(1,i), GRAVITY_PARAMS(2,i), GRAVITY_PARAMS(3,i)]
                CALL miyamotonagaipotential(p3, N, x, y, z, phi_c)
            CASE (GRAVITY_KIND_LONGMURALIBAR)
                p4 = [GRAVITY_PARAMS(1,i), GRAVITY_PARAMS(2,i), &
                      GRAVITY_PARAMS(3,i), GRAVITY_PARAMS(4,i)]
                CALL longmuralibarpotential(p4, N, x, y, z, phi_c)
            CASE DEFAULT
                WRITE(*,'(A)') "WARNING: evaluategravitypotential: unknown component kind"
            END SELECT
            phi = phi + phi_c
        END DO
    END SUBROUTINE evaluategravitypotential

    SUBROUTINE printgravitystate()
        IMPLICIT NONE
        INTEGER :: i
        CHARACTER(LEN=32) :: kind_name

        WRITE(*,'(A)') "=== gravity module state ==="
        WRITE(*,'(A,L1)')    "  GRAVITY_FINALIZED    : ", GRAVITY_FINALIZED
        WRITE(*,'(A,I0)')    "  GRAVITY_NCOMP        : ", GRAVITY_NCOMP
        IF (GRAVITY_G_IS_DEFAULT) THEN
            WRITE(*,'(A,ES14.8,A)') "  GRAVITY_G            : ", GRAVITY_G, " (default)"
        ELSE
            WRITE(*,'(A,ES14.8,A)') "  GRAVITY_G            : ", GRAVITY_G, " (user override)"
        END IF
        DO i = 1, GRAVITY_NCOMP
            SELECT CASE (GRAVITY_KIND(i))
            CASE (GRAVITY_KIND_PLUMMER);        kind_name = "plummer"
            CASE (GRAVITY_KIND_HERNQUIST);      kind_name = "hernquist"
            CASE (GRAVITY_KIND_MIYAMOTONAGAI);  kind_name = "miyamotonagai"
            CASE (GRAVITY_KIND_LONGMURALIBAR);  kind_name = "longmuralibar"
            CASE DEFAULT;                        kind_name = "unknown"
            END SELECT
            WRITE(*,'(A,I0,A,A)') "  component ", i, ": ", TRIM(kind_name)
        END DO
        WRITE(*,'(A)') "==========================="
    END SUBROUTINE printgravitystate

    
    ! =======================================================================
    ! =======================================================================
    ! =======================================================================
    ! =======================================================================
    ! ANALYTICAL POTENTIALS
    ! Closed-form gravity models with direct evaluation
    ! =======================================================================
    ! =======================================================================
    ! =======================================================================
    ! =======================================================================


    SUBROUTINE hernquistforce(params, N, x, y, z, force)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(2) :: params
        REAL*8, INTENT(OUT), DIMENSION(N,3) :: force
        REAL*8, DIMENSION(N) :: r, amod
        REAL*8 :: M, a

        M = params(1)
        a = params(2)
        r = SQRT(x*x + y*y + z*z)
        amod = -GRAVITY_G*M / (r + a)**2

        force(:,1) = amod*x
        force(:,2) = amod*y
        force(:,3) = amod*z
    END SUBROUTINE hernquistforce

    SUBROUTINE hernquistpotential(params, N, x, y, z, phi)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(2) :: params
        REAL*8, INTENT(OUT), DIMENSION(N) :: phi
        REAL*8, DIMENSION(N) :: r
        REAL*8 :: M, a

        M = params(1)
        a = params(2)
        r = SQRT(x*x + y*y + z*z)
        phi = -GRAVITY_G*M / (r + a)
    END SUBROUTINE hernquistpotential

    SUBROUTINE plummerforce(params, N, x, y, z, force)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(2) :: params
        REAL*8, INTENT(OUT), DIMENSION(N,3) :: force
        REAL*8, DIMENSION(N) :: r, amod
        REAL*8 :: M, b

        M = params(1)
        b = params(2)
        r = SQRT(x*x + y*y + z*z)
        amod = -GRAVITY_G*M / (r*r + b*b)**1.5

        force(:,1) = amod*x
        force(:,2) = amod*y
        force(:,3) = amod*z
    END SUBROUTINE plummerforce

    SUBROUTINE plummerpotential(params, N, x, y, z, phi)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(2) :: params
        REAL*8, INTENT(OUT), DIMENSION(N) :: phi
        REAL*8, DIMENSION(N) :: r
        REAL*8 :: M, b

        M = params(1)
        b = params(2)
        r = SQRT(x*x + y*y + z*z)
        phi = -GRAVITY_G*M / SQRT(r*r + b*b)
    END SUBROUTINE plummerpotential

    SUBROUTINE longmuralibarforce(params, N, x, y, z, force)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(4) :: params
        REAL*8, INTENT(OUT), DIMENSION(N,3) :: force
        REAL*8 :: M, abar, bbar, cbar
        REAL*8, DIMENSION(N) :: Tplus, Tminus

        M = params(1)
        abar = params(2)
        bbar = params(3)
        cbar = params(4)

        Tplus = SQRT((abar+x)**2 + y*y + (bbar+SQRT(cbar*cbar+z*z))**2)
        Tminus = SQRT((abar-x)**2 + y*y + (bbar+SQRT(cbar*cbar+z*z))**2)

        force(:,1) = -2.0D0*GRAVITY_G*M*x / ((Tplus*Tminus)*(Tplus+Tminus))
        force(:,2) = -GRAVITY_G*M*y/((2.0D0*Tplus*Tminus)*(y*y+(bbar+SQRT(z*z+cbar*cbar))**2)) * &
                     (Tplus+Tminus-4.0D0*x*x/(Tplus+Tminus))
        force(:,3) = -GRAVITY_G*M*z/((2.0D0*Tplus*Tminus)*(y*y+(bbar+SQRT(z*z+cbar*cbar))**2)) * &
                     (Tplus+Tminus-4.0D0*x*x/(Tplus+Tminus)) * &
                     ((bbar+SQRT(z*z+cbar*cbar))/SQRT(z*z+cbar*cbar))
    END SUBROUTINE longmuralibarforce

    SUBROUTINE longmuralibarpotential(params, N, x, y, z, phi)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(4) :: params
        REAL*8, INTENT(OUT), DIMENSION(N) :: phi
        REAL*8 :: M, abar, bbar, cbar
        REAL*8, DIMENSION(N) :: Tplus, Tminus

        M = params(1)
        abar = params(2)
        bbar = params(3)
        cbar = params(4)

        Tplus = SQRT((abar+x)**2 + y*y + (bbar+SQRT(cbar*cbar+z*z))**2)
        Tminus = SQRT((abar-x)**2 + y*y + (bbar+SQRT(cbar*cbar+z*z))**2)
        phi = (GRAVITY_G*M / (2.0D0*abar)) * LOG((x-abar+Tminus)/(x+abar+Tplus))
    END SUBROUTINE longmuralibarpotential

    SUBROUTINE allensantillianhaloforce(params, N, x, y, z, force)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(4) :: params
        REAL*8, INTENT(OUT), DIMENSION(N,3) :: force
        REAL*8 :: M, scale_length, exp, cutoffradius, Mtot, dcut
        REAL*8, DIMENSION(N) :: r, amod, d, d_exp_minus_1, d_exp_minus_3
        LOGICAL, DIMENSION(N) :: outside_cutoff, at_zero

        M = params(1)
        scale_length = params(2)
        exp = params(3)
        cutoffradius = params(4)

        r = SQRT(x*x + y*y + z*z)
        d = r / scale_length
        dcut = cutoffradius / scale_length
        d_exp_minus_1 = d**(exp-1)
        d_exp_minus_3 = d**(exp-3)

        Mtot = M * dcut**exp / (1.0D0 + dcut**(exp-1))
        amod = -(GRAVITY_G*M/scale_length**3) * (d_exp_minus_3 / (1.0D0 + d_exp_minus_1))

        outside_cutoff = (r > cutoffradius)
        at_zero = (r == 0.0D0)
        amod = MERGE(-(GRAVITY_G*Mtot/r**3), amod, outside_cutoff)
        amod = MERGE(0.0D0, amod, at_zero)

        force(:,1) = amod*x
        force(:,2) = amod*y
        force(:,3) = amod*z
    END SUBROUTINE allensantillianhaloforce

    SUBROUTINE allensantillianhalopotential(params, N, x, y, z, phi)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(4) :: params
        REAL*8, INTENT(OUT), DIMENSION(N) :: phi
        REAL*8, DIMENSION(N) :: term1, r, d, d_exp_minus_1
        REAL*8 :: M, scale_length, exp, cutoffradius, Mtot, term2, dcut
        LOGICAL, DIMENSION(N) :: outside_cutoff

        M = params(1)
        scale_length = params(2)
        exp = params(3)
        cutoffradius = params(4)

        r = SQRT(x*x + y*y + z*z)
        d = r / scale_length
        dcut = cutoffradius / scale_length
        d_exp_minus_1 = d**(exp-1)
        term1 = 1.0D0 + d_exp_minus_1
        term2 = 1.0D0 + dcut**(exp-1)

        Mtot = M * dcut**exp / (1.0D0 + dcut**(exp-1))
        phi = GRAVITY_G*M/(scale_length*(exp-1.0D0)) * LOG(term1/term2) - GRAVITY_G*Mtot/cutoffradius

        outside_cutoff = (r > cutoffradius)
        phi = MERGE(-GRAVITY_G*Mtot/r, phi, outside_cutoff)
    END SUBROUTINE allensantillianhalopotential

    SUBROUTINE miyamotonagaiforce(params, N, x, y, z, force)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(3) :: params
        REAL*8, INTENT(OUT), DIMENSION(N,3) :: force
        REAL*8, DIMENSION(N) :: R, amod, zmod
        REAL*8 :: M, a, b

        M = params(1)
        a = params(2)
        b = params(3)

        R = SQRT(x*x + y*y)
        zmod = a + SQRT(z*z + b*b)
        amod = -GRAVITY_G*M / (R*R + zmod*zmod)**1.5

        force(:,1) = amod*x
        force(:,2) = amod*y
        force(:,3) = amod*z*zmod/SQRT(z*z+b*b)
    END SUBROUTINE miyamotonagaiforce

    SUBROUTINE miyamotonagaipotential(params, N, x, y, z, phi)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(3) :: params
        REAL*8, INTENT(OUT), DIMENSION(N) :: phi
        REAL*8, DIMENSION(N) :: R, zmod
        REAL*8 :: M, a, b

        M = params(1)
        a = params(2)
        b = params(3)
        R = SQRT(x*x + y*y)
        zmod = a + SQRT(z*z + b*b)
        phi = -GRAVITY_G*M / SQRT(R*R + zmod*zmod)
    END SUBROUTINE miyamotonagaipotential

    SUBROUTINE pouliasis2017piiforce(params, N, x, y, z, force)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(10) :: params
        REAL*8, INTENT(OUT), DIMENSION(N,3) :: force
        REAL*8, DIMENSION(3) :: thindisk, thickdisk
        REAL*8, DIMENSION(4) :: halo
        REAL*8, DIMENSION(N,3) :: force_h, force_d1, force_d2

        halo = (/params(1), params(2), params(3), params(4)/)
        thindisk = (/params(5), params(6), params(7)/)
        thickdisk = (/params(8), params(9), params(10)/)

        CALL allensantillianhaloforce(halo, N, x, y, z, force_h)
        CALL miyamotonagaiforce(thindisk, N, x, y, z, force_d1)
        CALL miyamotonagaiforce(thickdisk, N, x, y, z, force_d2)
        force = force_h + force_d1 + force_d2
    END SUBROUTINE pouliasis2017piiforce

    SUBROUTINE pouliasis2017piipotential(params, N, x, y, z, phi)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(10) :: params
        REAL*8, INTENT(OUT), DIMENSION(N) :: phi
        REAL*8, DIMENSION(3) :: thindisk, thickdisk
        REAL*8, DIMENSION(4) :: halo
        REAL*8, DIMENSION(N) :: phi_h, phi_d1, phi_d2

        halo = (/params(1), params(2), params(3), params(4)/)
        thindisk = (/params(5), params(6), params(7)/)
        thickdisk = (/params(8), params(9), params(10)/)

        CALL allensantillianhalopotential(halo, N, x, y, z, phi_h)
        CALL miyamotonagaipotential(thindisk, N, x, y, z, phi_d1)
        CALL miyamotonagaipotential(thickdisk, N, x, y, z, phi_d2)
        phi = phi_h + phi_d1 + phi_d2
    END SUBROUTINE pouliasis2017piipotential

    SUBROUTINE NBODYPLUMMERS(params,N,x,y,z,ax,ay,az,phiTensor)
        ! Computeres the inter gravitational forces between N particles
        ! uses a plummer sphere for each particle. i.e. individual softening parameters for each body... 
        IMPLICIT NONE 
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN),DIMENSION(2*N + 1) :: params  ! first is the gravitational constant, then the masses, then the radii
        REAL*8, INTENT(IN),DIMENSION(N) :: x,y,z
        REAL*8, INTENT(OUT),DIMENSION(N) :: ax,ay,az
        REAL*8, INTENT(OUT),DIMENSION(N,N) :: phiTensor
        REAL*8, DIMENSION(N,N) :: FX,FY,FZ ! the forces on each particle
        REAL*8, DIMENSION(N) :: masses,scaleradii 
        REAL*8, DIMENSION(2) :: params2
        REAL*8, DIMENSION(N,3) :: force
        integer :: i,j
        masses=params(1+1:N+1)
        scaleradii=params(N+1+1:2*N+1) ! plus ones for the gravitational constant 
        ! force from i on all the others
        DO j=1,N
            params2(1)=masses(j)
            params2(2)=scaleradii(j)
            ! potential calculates the force of i on all the other particles
            CALL plummerforce(params2, N, x-x(j), y-y(j), z-z(j), force)
            CALL plummerpotential(params2, N, x-x(j), y-y(j), z-z(j), phiTensor(:,j))
            ax = force(:,1)
            ay = force(:,2)
            az = force(:,3)
            FX(:,j)=ax*masses
            FY(:,j)=ay*masses
            FZ(:,j)=az*masses
            phiTensor(j,j)=0
            FX(j,j)=0
            FY(j,j)=0
            FZ(j,j)=0
        END DO
        ! now sum up the forces
        ! we don't add a minus sign because we sum along the rows, which means force ON i from the others... 
        DO i=1,N
            ax(i)=sum(FX(i,:))/masses(i)
            ay(i)=sum(FY(i,:))/masses(i)
            az(i)=sum(FZ(i,:))/masses(i)
        END DO
    END SUBROUTINE NBODYPLUMMERS   
    
    SUBROUTINE pointmassconfiguration(NParticles,masses,xGC,yGC,zGC,N,x,y,z,ax,ay,az,phi)
        ! given a configuration of point masses,
        ! find the acceleration and potential at a given point
        IMPLICIT NONE 
        INTEGER, INTENT(IN) :: N,NParticles
        REAL*8, INTENT(IN),DIMENSION(N) :: x,y,z
        REAL*8, INTENT(OUT),DIMENSION(N) :: ax,ay,az
        REAL*8, INTENT(OUT),DIMENSION(N) :: phi
        REAL*8, INTENT(IN),DIMENSION(NParticles) :: xGC, yGC, zGC, masses
        REAL*8 :: dx,dy,dz,dr,dr3
        INTEGER :: i,j

        ! initialize at zero
        ax=0
        ay=0
        az=0
        DO i=1,N
            DO j=1,NParticles
                dx=xGC(j)-x(i)
                dy=yGC(j)-y(i)
                dz=zGC(j)-z(i)
                dr=sqrt(dx*dx+dy*dy+dz*dz)
                dr3=dr*dr*dr
                ax(i)=ax(i)+GRAVITY_G*masses(j)*dx/dr3
                ay(i)=ay(i)+GRAVITY_G*masses(j)*dy/dr3
                az(i)=az(i)+GRAVITY_G*masses(j)*dz/dr3
                phi(i)=phi(i)-GRAVITY_G*masses(j)/dr
            END DO
        END DO
    end SUBROUTINE pointmassconfiguration


    ! =======================================================================
    ! =======================================================================
    ! =======================================================================
    ! SPHERICAL HARMONICS INFRASTRUCTURE
    ! Legendre polynomial basis, projection, and evaluation primitives
    ! =======================================================================
    ! =======================================================================
    ! =======================================================================

    SUBROUTINE initaxisymmetricbasisexpansion(lmax, nr, r_grid)
        ! Allocate basis-expansion storage and store the radial grid.
        ! Does NOT project any density or compute potential tables.
        ! Call a density subroutine (e.g. exponential_oblate_halo) afterwards
        ! to trigger projection and table construction.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: lmax, nr
        REAL*8, INTENT(IN), DIMENSION(nr) :: r_grid

        CALL clearaxisymmetricbasisexpansion()

        ALLOCATE(BASIS_R_GRID(nr))
        ALLOCATE(BASIS_RHO_L_GRID(0:lmax, nr))
        ALLOCATE(BASIS_PHI_L_GRID(0:lmax, nr))
        ALLOCATE(BASIS_DPHI_L_DR_GRID(0:lmax, nr))

        BASIS_RHO_L_GRID     = 0.0D0
        BASIS_PHI_L_GRID     = 0.0D0
        BASIS_DPHI_L_DR_GRID = 0.0D0

        BASIS_G    = GRAVITY_G
        BASIS_LMAX = lmax
        BASIS_NR   = nr
        BASIS_R_GRID = r_grid

        BASIS_GRID_SET              = .TRUE.
        BASIS_EXPANSION_INITIALIZED = .FALSE.
    END SUBROUTINE initaxisymmetricbasisexpansion

    SUBROUTINE clearaxisymmetricbasisexpansion()
        IMPLICIT NONE

        IF (ALLOCATED(BASIS_R_GRID))          DEALLOCATE(BASIS_R_GRID)
        IF (ALLOCATED(BASIS_RHO_L_GRID))      DEALLOCATE(BASIS_RHO_L_GRID)
        IF (ALLOCATED(BASIS_PHI_L_GRID))      DEALLOCATE(BASIS_PHI_L_GRID)
        IF (ALLOCATED(BASIS_DPHI_L_DR_GRID))  DEALLOCATE(BASIS_DPHI_L_DR_GRID)

        BASIS_GRID_SET              = .FALSE.
        BASIS_EXPANSION_INITIALIZED = .FALSE.
        BASIS_LMAX = -1
        BASIS_NR   = -1
        BASIS_G    = -1.0D0
    END SUBROUTINE clearaxisymmetricbasisexpansion

    SUBROUTINE default_init_basis_expansion()
        ! Auto-initialize with sensible defaults when the user has not called
        ! initaxisymmetricbasisexpansion explicitly.
        ! Grid: 100 log-spaced points from 1e-4 to 1e3; lmax = 20.
        IMPLICIT NONE
        INTEGER, PARAMETER :: default_lmax = 20
        INTEGER, PARAMETER :: default_nr   = 100
        REAL*8,  PARAMETER :: default_rmin = 1.0D-4
        REAL*8,  PARAMETER :: default_rmax = 1.0D3
        REAL*8 :: r_grid(default_nr)
        REAL*8 :: log_rmin, log_rmax, dlog_r
        INTEGER :: i

        log_rmin = LOG(default_rmin)
        log_rmax = LOG(default_rmax)
        dlog_r   = (log_rmax - log_rmin) / DBLE(default_nr - 1)
        DO i = 1, default_nr
            r_grid(i) = EXP(log_rmin + (i-1) * dlog_r)
        END DO
        CALL initaxisymmetricbasisexpansion(default_lmax, default_nr, r_grid)
    END SUBROUTINE default_init_basis_expansion    

    SUBROUTINE compute_phi_tables_from_rho()
        ! Build potential tables BASIS_PHI_L_GRID and BASIS_DPHI_L_DR_GRID from
        ! the projected density coefficients BASIS_RHO_L_GRID via the Green's
        ! function integrals for each l-mode:
        !   Phi_l(r) = -4*pi*G/(2l+1) * [ r^{-(l+1)} I_l^<(r) + r^l I_l^>(r) ]
        ! where I_l^<(r) = int_0^r  r'^{l+2}   rho_l(r') dr'
        !       I_l^>(r) = int_r^inf r'^{1-l} rho_l(r') dr'
        ! dPhi_l/dr = -4*pi*G/(2l+1) * [ -(l+1)*r^{-(l+2)} I_l^< + l*r^{l-1} I_l^> ]
        IMPLICIT NONE
        REAL*8, PARAMETER :: pi_phi = 3.14159265358979323846D0
        INTEGER :: l, i
        REAL*8  :: prefactor, r_lo, r_hi, f_lo, f_hi, dr
        REAL*8, ALLOCATABLE :: I_less(:), I_greater(:)

        ALLOCATE(I_less(BASIS_NR), I_greater(BASIS_NR))

        DO l = 0, BASIS_LMAX, 2
            prefactor = -4.0D0 * pi_phi * BASIS_G / DBLE(2*l + 1)

            ! Forward cumulative trapezoid: I_less(i) = int_0^r r'^{l+2} rho_l dr'
            I_less(1) = 0.0D0
            DO i = 2, BASIS_NR
                r_lo = BASIS_R_GRID(i-1);  r_hi = BASIS_R_GRID(i)
                f_lo = r_lo**(l+2) * BASIS_RHO_L_GRID(l, i-1)
                f_hi = r_hi**(l+2) * BASIS_RHO_L_GRID(l, i)
                I_less(i) = I_less(i-1) + 0.5D0 * (r_hi - r_lo) * (f_lo + f_hi)
            END DO

            ! Backward cumulative trapezoid: I_greater(i) = int_r^inf r'^{1-l} rho_l dr'
            I_greater(BASIS_NR) = 0.0D0
            DO i = BASIS_NR-1, 1, -1
                r_lo = BASIS_R_GRID(i);  r_hi = BASIS_R_GRID(i+1)
                f_lo = r_lo**(1-l) * BASIS_RHO_L_GRID(l, i)
                f_hi = r_hi**(1-l) * BASIS_RHO_L_GRID(l, i+1)
                dr = r_hi - r_lo
                I_greater(i) = I_greater(i+1) + 0.5D0 * dr * (f_lo + f_hi)
            END DO

            DO i = 1, BASIS_NR
                BASIS_PHI_L_GRID(l, i) = prefactor * ( &
                    BASIS_R_GRID(i)**(-(l+1)) * I_less(i) + &
                    BASIS_R_GRID(i)**l        * I_greater(i) )
                BASIS_DPHI_L_DR_GRID(l, i) = prefactor * ( &
                    -(l+1) * BASIS_R_GRID(i)**(-(l+2)) * I_less(i) + &
                    DBLE(l) * BASIS_R_GRID(i)**(l-1)   * I_greater(i) )
            END DO
        END DO

        DEALLOCATE(I_less, I_greater)
    END SUBROUTINE compute_phi_tables_from_rho

    SUBROUTINE axisymmetricbasisexpansion_eval(N, x, y, z, ax, ay, az, phi_out)
        ! Evaluate accelerations and potential for all N particles by
        ! interpolating the pre-computed l-mode tables and summing over modes.
        ! Cartesian force: a_i = -d_i phi, chain rule via (r, mu=z/r).
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN),  DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi_out
        INTEGER :: i, l, jlo, jhi, jmid
        REAL*8  :: r, r_safe, mu, R_cyl_sq
        REAL*8  :: alpha_r, phi_l_r, dphi_l_dr_r
        REAL*8  :: phi_v, dphi_dr_v, dphi_dmu_v
        REAL*8  :: p(0:BASIS_LMAX), dp_dmu(0:BASIS_LMAX)
        REAL*8, PARAMETER :: eps_r = 1.0D-30

        DO i = 1, N
            r      = SQRT(x(i)**2 + y(i)**2 + z(i)**2)
            r_safe = MAX(r, eps_r)
            mu     = z(i) / r_safe
            R_cyl_sq = x(i)**2 + y(i)**2

            ! Binary search for bracketing index in r_grid
            IF (r_safe <= BASIS_R_GRID(1)) THEN
                jlo = 1;  alpha_r = 0.0D0
            ELSE IF (r_safe >= BASIS_R_GRID(BASIS_NR)) THEN
                jlo = BASIS_NR - 1;  alpha_r = 1.0D0
            ELSE
                jlo = 1;  jhi = BASIS_NR
                DO WHILE (jhi - jlo > 1)
                    jmid = (jlo + jhi) / 2
                    IF (BASIS_R_GRID(jmid) <= r_safe) THEN
                        jlo = jmid
                    ELSE
                        jhi = jmid
                    END IF
                END DO
                ! Log-r interpolation weight (accurate on log-spaced grids)
                alpha_r = LOG(r_safe / BASIS_R_GRID(jlo)) / &
                          LOG(BASIS_R_GRID(jlo+1) / BASIS_R_GRID(jlo))
            END IF

            CALL legendre_axisymmetric_basis(BASIS_LMAX, mu, p, dp_dmu)

            phi_v = 0.0D0;  dphi_dr_v = 0.0D0;  dphi_dmu_v = 0.0D0
            DO l = 0, BASIS_LMAX, 2
                phi_l_r     = linear_interp_scalar(BASIS_PHI_L_GRID(l,jlo),     BASIS_PHI_L_GRID(l,jlo+1),     alpha_r)
                dphi_l_dr_r = linear_interp_scalar(BASIS_DPHI_L_DR_GRID(l,jlo), BASIS_DPHI_L_DR_GRID(l,jlo+1), alpha_r)
                phi_v      = phi_v      + phi_l_r     * p(l)
                dphi_dr_v  = dphi_dr_v  + dphi_l_dr_r * p(l)
                dphi_dmu_v = dphi_dmu_v + phi_l_r     * dp_dmu(l)
            END DO

            phi_out(i) = phi_v
            ! a = -grad(phi).  With phi(r,mu), mu=z/r:
            !   d_phi/dx = phi_r*(x/r) - phi_mu*(z*x/r^3)
            !   d_phi/dy = phi_r*(y/r) - phi_mu*(z*y/r^3)
            !   d_phi/dz = phi_r*(z/r) + phi_mu*(R^2/r^3)
            ax(i) = -dphi_dr_v * x(i)/r_safe + dphi_dmu_v * z(i)*x(i)/r_safe**3
            ay(i) = -dphi_dr_v * y(i)/r_safe + dphi_dmu_v * z(i)*y(i)/r_safe**3
            az(i) = -dphi_dr_v * z(i)/r_safe - dphi_dmu_v * R_cyl_sq/r_safe**3
        END DO
    END SUBROUTINE axisymmetricbasisexpansion_eval

    SUBROUTINE compute_phi_tables_from_rho_component(r_grid, rho_l_grid, phi_l_grid, dphi_l_dr_grid)
        IMPLICIT NONE
        REAL*8, INTENT(IN), DIMENSION(:) :: r_grid
        REAL*8, INTENT(IN), DIMENSION(:,:) :: rho_l_grid
        REAL*8, INTENT(OUT), DIMENSION(:,:) :: phi_l_grid, dphi_l_dr_grid
        REAL*8, PARAMETER :: pi_phi = 3.14159265358979323846D0
        INTEGER :: l, i, l_lo, l_hi, nr
        REAL*8  :: prefactor, r_lo, r_hi, f_lo, f_hi, dr
        REAL*8, ALLOCATABLE :: I_less(:), I_greater(:)

        nr = SIZE(r_grid)
        l_lo = LBOUND(rho_l_grid, 1)
        l_hi = UBOUND(rho_l_grid, 1)

        ALLOCATE(I_less(nr), I_greater(nr))

        DO l = l_lo, l_hi, 2
            prefactor = -4.0D0 * pi_phi * GRAVITY_G / DBLE(2*l + 1)

            I_less(1) = 0.0D0
            DO i = 2, nr
                r_lo = r_grid(i-1);  r_hi = r_grid(i)
                f_lo = r_lo**(l+2) * rho_l_grid(l, i-1)
                f_hi = r_hi**(l+2) * rho_l_grid(l, i)
                I_less(i) = I_less(i-1) + 0.5D0 * (r_hi - r_lo) * (f_lo + f_hi)
            END DO

            I_greater(nr) = 0.0D0
            DO i = nr-1, 1, -1
                r_lo = r_grid(i);  r_hi = r_grid(i+1)
                f_lo = r_lo**(1-l) * rho_l_grid(l, i)
                f_hi = r_hi**(1-l) * rho_l_grid(l, i+1)
                dr = r_hi - r_lo
                I_greater(i) = I_greater(i+1) + 0.5D0 * dr * (f_lo + f_hi)
            END DO

            DO i = 1, nr
                phi_l_grid(l, i) = prefactor * ( &
                    r_grid(i)**(-(l+1)) * I_less(i) + &
                    r_grid(i)**l        * I_greater(i) )
                dphi_l_dr_grid(l, i) = prefactor * ( &
                    -(l+1) * r_grid(i)**(-(l+2)) * I_less(i) + &
                    DBLE(l) * r_grid(i)**(l-1)   * I_greater(i) )
            END DO
        END DO

        DEALLOCATE(I_less, I_greater)
    END SUBROUTINE compute_phi_tables_from_rho_component

    SUBROUTINE axisymmetricbasisexpansion_eval_component(N, x, y, z, r_grid, phi_l_grid, dphi_l_dr_grid, ax, ay, az, phi_out)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN),  DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN),  DIMENSION(:) :: r_grid
        REAL*8, INTENT(IN),  DIMENSION(:,:) :: phi_l_grid, dphi_l_dr_grid
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi_out
        INTEGER :: i, l, jlo, jhi, jmid, l_lo, l_hi, nr
        REAL*8  :: r, r_safe, mu, R_cyl_sq
        REAL*8  :: alpha_r, phi_l_r, dphi_l_dr_r
        REAL*8  :: phi_v, dphi_dr_v, dphi_dmu_v
        REAL*8, ALLOCATABLE :: p(:), dp_dmu(:)
        REAL*8, PARAMETER :: eps_r = 1.0D-30

        nr = SIZE(r_grid)
        l_lo = LBOUND(phi_l_grid, 1)
        l_hi = UBOUND(phi_l_grid, 1)
        ALLOCATE(p(l_lo:l_hi), dp_dmu(l_lo:l_hi))

        DO i = 1, N
            r      = SQRT(x(i)**2 + y(i)**2 + z(i)**2)
            r_safe = MAX(r, eps_r)
            mu     = z(i) / r_safe
            R_cyl_sq = x(i)**2 + y(i)**2

            IF (r_safe <= r_grid(1)) THEN
                jlo = 1;  alpha_r = 0.0D0
            ELSE IF (r_safe >= r_grid(nr)) THEN
                jlo = nr - 1;  alpha_r = 1.0D0
            ELSE
                jlo = 1;  jhi = nr
                DO WHILE (jhi - jlo > 1)
                    jmid = (jlo + jhi) / 2
                    IF (r_grid(jmid) <= r_safe) THEN
                        jlo = jmid
                    ELSE
                        jhi = jmid
                    END IF
                END DO
                alpha_r = LOG(r_safe / r_grid(jlo)) / LOG(r_grid(jlo+1) / r_grid(jlo))
            END IF

            CALL legendre_axisymmetric_basis(l_hi, mu, p, dp_dmu)

            phi_v = 0.0D0;  dphi_dr_v = 0.0D0;  dphi_dmu_v = 0.0D0
            DO l = l_lo, l_hi, 2
                phi_l_r     = linear_interp_scalar(phi_l_grid(l,jlo),     phi_l_grid(l,jlo+1),     alpha_r)
                dphi_l_dr_r = linear_interp_scalar(dphi_l_dr_grid(l,jlo), dphi_l_dr_grid(l,jlo+1), alpha_r)
                phi_v      = phi_v      + phi_l_r     * p(l)
                dphi_dr_v  = dphi_dr_v  + dphi_l_dr_r * p(l)
                dphi_dmu_v = dphi_dmu_v + phi_l_r     * dp_dmu(l)
            END DO

            phi_out(i) = phi_v
            ax(i) = -dphi_dr_v * x(i)/r_safe + dphi_dmu_v * z(i)*x(i)/r_safe**3
            ay(i) = -dphi_dr_v * y(i)/r_safe + dphi_dmu_v * z(i)*y(i)/r_safe**3
            az(i) = -dphi_dr_v * z(i)/r_safe - dphi_dmu_v * R_cyl_sq/r_safe**3
        END DO

        DEALLOCATE(p, dp_dmu)
    END SUBROUTINE axisymmetricbasisexpansion_eval_component

    ! =======================================================================
    ! =======================================================================
    ! =======================================================================
    ! BESSEL/TABLE INFRASTRUCTURE
    ! Disk-table construction via Bessel/Hankel quadrature and table evaluation
    ! =======================================================================
    ! =======================================================================
    ! =======================================================================

    SUBROUTINE exponential_disk_bessel_eval_component(params, N, x, y, z, ax, ay, az, phi)
        ! Exponential-disk application of the generic Bessel/Poisson method.
        ! This is the preferred basis path for strongly flattened disk-like
        ! components (e.g. q <= 0.3), while Legendre is better for near-
        ! spherical components.
        !
        ! Density model associated with this evaluator:
        !   rho(R,z) = sigma0/(2*hZ) * exp(-R/hR - |z|/hZ)
        IMPLICIT NONE
        REAL*8, INTENT(IN), DIMENSION(*) :: params
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi
        INTEGER, PARAMETER :: nk = 256
        REAL*8, PARAMETER :: pi = 3.14159265358979323846D0
        REAL*8, PARAMETER :: eps = 1.0D-16
        INTEGER :: i, j
        REAL*8 :: R, absz, signz, aR, A
        REAL*8 :: sigma0, hR, hZ
        REAL*8 :: sum_phi, sum_ar, sum_az
        REAL*8 :: kernel, j0, j1
        REAL*8 :: kval, t, mu
        REAL*8, DIMENSION(nk) :: k_grid, quad_w, base_kernel, mu_q, w_q

        sigma0 = params(1)
        hR = params(2)
        hZ = params(3)

        A = 2.0D0 * pi * GRAVITY_G * sigma0 * hR * hR

        ! Integrate k in [0, +inf) via t in [0,1):
        !   k = t/(1-t),   dk = dt/(1-t)^2
        ! This avoids a hard k cutoff and substantially reduces ringing.
        CALL gauss_legendre_nodes_weights(nk, mu_q, w_q)
        DO j = 1, nk
            mu = mu_q(j)
            t = 0.5D0 * (mu + 1.0D0)
            kval = t / MAX(1.0D0 - t, 1.0D-14)
            k_grid(j) = kval

            ! 0.5 converts [-1,1] -> [0,1], jacobian handles [0,1) -> [0,+inf)
            quad_w(j) = 0.5D0 * w_q(j) / MAX((1.0D0 - t)**2, 1.0D-14)
            base_kernel(j) = 1.0D0 / ( (1.0D0 + (kval*hR)**2)**1.5D0 * (1.0D0 + kval*hZ) )
        END DO

        DO i = 1, N
            R = SQRT(x(i)**2 + y(i)**2)
            absz = ABS(z(i))
            IF (z(i) > 0.0D0) THEN
                signz = 1.0D0
            ELSE IF (z(i) < 0.0D0) THEN
                signz = -1.0D0
            ELSE
                signz = 0.0D0
            END IF

            sum_phi = 0.0D0
            sum_ar = 0.0D0
            sum_az = 0.0D0
            DO j = 1, nk
                kval = k_grid(j)
                kernel = EXP(-kval * absz) * base_kernel(j)
                j0 = bessel_j0_scalar(kval * R)
                j1 = bessel_j1_scalar(kval * R)

                sum_phi = sum_phi + quad_w(j) * j0 * kernel
                sum_ar  = sum_ar  + quad_w(j) * kval * j1 * kernel
                sum_az  = sum_az  + quad_w(j) * kval * j0 * kernel
            END DO

            phi(i) = -A * sum_phi
            aR = -A * sum_ar
            az(i) = -A * signz * sum_az

            IF (R > eps) THEN
                ax(i) = aR * x(i) / R
                ay(i) = aR * y(i) / R
            ELSE
                ax(i) = 0.0D0
                ay(i) = 0.0D0
            END IF
        END DO
    END SUBROUTINE exponential_disk_bessel_eval_component

    SUBROUTINE build_exponential_disk_table(component_index, sigma0, hR, hZ)
        ! Build a precomputed 2D (R, z>=0) potential table and derivative tables
        ! for conservative bicubic interpolation:
        !   Phi, dPhi/dR, dPhi/dz, d2Phi/(dR dz)
        ! Runtime forces are then computed as gradients of ONE interpolated Phi.
        !
        ! Table R extent: [0.001*hR, 200*hR], COMPOSITE_DISK_TABLE_NR nodes, log-spaced.
        ! Table z extent: [0, 200*hZ],        COMPOSITE_DISK_TABLE_NZ nodes, linear.
        ! Quadrature nodes: nk=512 (higher accuracy acceptable for offline build).
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: component_index
        REAL*8, INTENT(IN) :: sigma0, hR, hZ
        INTEGER, PARAMETER :: nk = 512
        REAL*8, PARAMETER :: pi = 3.14159265358979323846D0
        INTEGER :: iR, iZ, j
        REAL*8 :: R, z_val, A
        REAL*8 :: logR_min, logR_max, dlogR, z_max, dz
        REAL*8 :: sum_phi, sum_dphi_dR, sum_dphi_dz, kernel, j0v, j1v
        REAL*8 :: kval, t, mu
        REAL*8, DIMENSION(nk) :: k_grid, quad_w, base_kernel, mu_q, w_q

        A = 2.0D0 * pi * GRAVITY_G * sigma0 * hR * hR

        logR_min = LOG(1.0D-3 * hR)
        logR_max = LOG(2.0D2 * hR)
        dlogR    = (logR_max - logR_min) / DBLE(COMPOSITE_DISK_TABLE_NR - 1)
        z_max    = 2.0D2 * hZ
        dz       = z_max / DBLE(COMPOSITE_DISK_TABLE_NZ - 1)

        ! Store grid metadata for runtime interpolation
        COMPOSITE_DISK_TABLE_META(1, component_index) = logR_min
        COMPOSITE_DISK_TABLE_META(2, component_index) = dlogR
        COMPOSITE_DISK_TABLE_META(3, component_index) = dz

        ! Precompute k quadrature nodes and k-independent kernel factor
        CALL gauss_legendre_nodes_weights(nk, mu_q, w_q)
        DO j = 1, nk
            mu        = mu_q(j)
            t         = 0.5D0 * (mu + 1.0D0)
            kval      = t / MAX(1.0D0 - t, 1.0D-14)
            k_grid(j) = kval
            quad_w(j) = 0.5D0 * w_q(j) / MAX((1.0D0 - t)**2, 1.0D-14)
            base_kernel(j) = 1.0D0 / ((1.0D0 + (kval*hR)**2)**1.5D0 * (1.0D0 + kval*hZ))
        END DO

        ! Evaluate Phi and first derivatives on each (R, z>=0) grid point.
        DO iR = 1, COMPOSITE_DISK_TABLE_NR
            R = EXP(logR_min + DBLE(iR - 1) * dlogR)
            DO iZ = 1, COMPOSITE_DISK_TABLE_NZ
                z_val = DBLE(iZ - 1) * dz

                sum_phi     = 0.0D0
                sum_dphi_dR = 0.0D0
                sum_dphi_dz = 0.0D0
                DO j = 1, nk
                    kval   = k_grid(j)
                    kernel = EXP(-kval * z_val) * base_kernel(j)
                    j0v    = bessel_j0_scalar(kval * R)
                    j1v    = bessel_j1_scalar(kval * R)

                    sum_phi     = sum_phi     + quad_w(j) * j0v * kernel
                    sum_dphi_dR = sum_dphi_dR + quad_w(j) * kval * j1v * kernel
                    sum_dphi_dz = sum_dphi_dz + quad_w(j) * kval * j0v * kernel
                END DO

                COMPOSITE_DISK_TABLE_PHI(iR, iZ, component_index)     = -A * sum_phi
                COMPOSITE_DISK_TABLE_DPHI_DR(iR, iZ, component_index) =  A * sum_dphi_dR
                COMPOSITE_DISK_TABLE_DPHI_DZ(iR, iZ, component_index) =  A * sum_dphi_dz
            END DO
        END DO

        ! Mixed derivative d/dz(dPhi/dR) from the precomputed first derivative.
        DO iR = 1, COMPOSITE_DISK_TABLE_NR
            DO iZ = 1, COMPOSITE_DISK_TABLE_NZ
                IF (iZ == 1) THEN
                    COMPOSITE_DISK_TABLE_D2PHI_DRDZ(iR, iZ, component_index) = &
                        (COMPOSITE_DISK_TABLE_DPHI_DR(iR, iZ+1, component_index) - &
                         COMPOSITE_DISK_TABLE_DPHI_DR(iR, iZ,   component_index)) / dz
                ELSE IF (iZ == COMPOSITE_DISK_TABLE_NZ) THEN
                    COMPOSITE_DISK_TABLE_D2PHI_DRDZ(iR, iZ, component_index) = &
                        (COMPOSITE_DISK_TABLE_DPHI_DR(iR, iZ,   component_index) - &
                         COMPOSITE_DISK_TABLE_DPHI_DR(iR, iZ-1, component_index)) / dz
                ELSE
                    COMPOSITE_DISK_TABLE_D2PHI_DRDZ(iR, iZ, component_index) = &
                        (COMPOSITE_DISK_TABLE_DPHI_DR(iR, iZ+1, component_index) - &
                         COMPOSITE_DISK_TABLE_DPHI_DR(iR, iZ-1, component_index)) / (2.0D0 * dz)
                END IF
            END DO
        END DO
    END SUBROUTINE build_exponential_disk_table

    SUBROUTINE disk_table_eval_component(component_index, N, x, y, z, ax, ay, az, phi)
        ! Conservative bicubic evaluation from one interpolated potential patch.
        ! Forces are computed as gradients of the same interpolated Phi:
        !   aR = -dPhi/dR, az = -dPhi/dz.
        ! Table stores z>=0 values; odd/even symmetry is applied for z<0.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: component_index, N
        REAL*8, INTENT(IN),  DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi
        REAL*8, PARAMETER :: eps = 1.0D-16
        INTEGER :: i, iR, iZ
        REAL*8 :: R, Rq, absz, signz, aR_val
        REAL*8 :: logR, logR_min, dlogR, dz, logR_max
        REAL*8 :: alphaR, alphaZ, zq
        REAL*8 :: R_lo, R_hi, dR, z_lo
        REAL*8 :: phi_loc, dphi_dR_loc, dphi_dz_abs
        REAL*8 :: f00, f10, f01, f11
        REAL*8 :: fx00, fx10, fx01, fx11
        REAL*8 :: fy00, fy10, fy01, fy11
        REAL*8 :: fxy00, fxy10, fxy01, fxy11

        logR_min = COMPOSITE_DISK_TABLE_META(1, component_index)
        dlogR    = COMPOSITE_DISK_TABLE_META(2, component_index)
        dz       = COMPOSITE_DISK_TABLE_META(3, component_index)
        logR_max = logR_min + DBLE(COMPOSITE_DISK_TABLE_NR - 1) * dlogR

        DO i = 1, N
            R    = SQRT(x(i)**2 + y(i)**2)
            absz = ABS(z(i))
            IF (z(i) > 0.0D0) THEN
                signz = 1.0D0
            ELSE IF (z(i) < 0.0D0) THEN
                signz = -1.0D0
            ELSE
                signz = 0.0D0
            END IF

            ! Clamp query coordinates to table domain.
            IF (R > eps) THEN
                logR = MIN(MAX(LOG(R), logR_min), logR_max)
                Rq = EXP(logR)
            ELSE
                logR = logR_min
                Rq = EXP(logR_min)
            END IF
            zq = MIN(absz, DBLE(COMPOSITE_DISK_TABLE_NZ - 1) * dz)

            ! Cell lookup in logR index space.
            iR = INT((logR - logR_min) / dlogR)
            iR = MAX(0, MIN(COMPOSITE_DISK_TABLE_NR - 2, iR))
            R_lo = EXP(logR_min + DBLE(iR) * dlogR)
            R_hi = EXP(logR_min + DBLE(iR + 1) * dlogR)
            dR = MAX(R_hi - R_lo, eps)
            alphaR = (Rq - R_lo) / dR
            alphaR = MAX(0.0D0, MIN(1.0D0, alphaR))

            ! Cell lookup in z index space.
            iZ = INT(zq / dz)
            iZ = MAX(0, MIN(COMPOSITE_DISK_TABLE_NZ - 2, iZ))
            z_lo = DBLE(iZ) * dz
            alphaZ = (zq - z_lo) / dz
            alphaZ = MAX(0.0D0, MIN(1.0D0, alphaZ))

            ! Gather patch values and derivatives at corners.
            f00   = COMPOSITE_DISK_TABLE_PHI(iR+1, iZ+1, component_index)
            f10   = COMPOSITE_DISK_TABLE_PHI(iR+2, iZ+1, component_index)
            f01   = COMPOSITE_DISK_TABLE_PHI(iR+1, iZ+2, component_index)
            f11   = COMPOSITE_DISK_TABLE_PHI(iR+2, iZ+2, component_index)

            fx00  = COMPOSITE_DISK_TABLE_DPHI_DR(iR+1, iZ+1, component_index)
            fx10  = COMPOSITE_DISK_TABLE_DPHI_DR(iR+2, iZ+1, component_index)
            fx01  = COMPOSITE_DISK_TABLE_DPHI_DR(iR+1, iZ+2, component_index)
            fx11  = COMPOSITE_DISK_TABLE_DPHI_DR(iR+2, iZ+2, component_index)

            fy00  = COMPOSITE_DISK_TABLE_DPHI_DZ(iR+1, iZ+1, component_index)
            fy10  = COMPOSITE_DISK_TABLE_DPHI_DZ(iR+2, iZ+1, component_index)
            fy01  = COMPOSITE_DISK_TABLE_DPHI_DZ(iR+1, iZ+2, component_index)
            fy11  = COMPOSITE_DISK_TABLE_DPHI_DZ(iR+2, iZ+2, component_index)

            fxy00 = COMPOSITE_DISK_TABLE_D2PHI_DRDZ(iR+1, iZ+1, component_index)
            fxy10 = COMPOSITE_DISK_TABLE_D2PHI_DRDZ(iR+2, iZ+1, component_index)
            fxy01 = COMPOSITE_DISK_TABLE_D2PHI_DRDZ(iR+1, iZ+2, component_index)
            fxy11 = COMPOSITE_DISK_TABLE_D2PHI_DRDZ(iR+2, iZ+2, component_index)

            CALL bicubic_hermite_eval_2d(f00, f10, f01, f11, &
                                         fx00, fx10, fx01, fx11, &
                                         fy00, fy10, fy01, fy11, &
                                         fxy00, fxy10, fxy01, fxy11, &
                                         dR, dz, alphaR, alphaZ, &
                                         phi_loc, dphi_dR_loc, dphi_dz_abs)

            phi(i) = phi_loc
            aR_val = -dphi_dR_loc
            az(i)  = -signz * dphi_dz_abs

            IF (R > eps) THEN
                ax(i) = aR_val * x(i) / R
                ay(i) = aR_val * y(i) / R
            ELSE
                ax(i) = 0.0D0
                ay(i) = 0.0D0
            END IF
        END DO
    END SUBROUTINE disk_table_eval_component

    ! =======================================================================
    ! =======================================================================
    ! =======================================================================
    ! DENSITY-ONLY POTENTIALS (Legendre-based)
    ! Models defined via density profiles, evaluated with Legendre expansion
    ! =======================================================================
    ! =======================================================================
    ! =======================================================================

    SUBROUTINE project_exponential_oblate_halo(rho0, s0, q)
        ! Project rho(r,mu) = rho0*exp(-r/s0*sqrt(1-(1-1/q^2)*mu^2)) onto
        ! Legendre modes using Gauss-Legendre quadrature in mu.
        ! Fills BASIS_RHO_L_GRID for even l only.
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: rho0, s0, q
        INTEGER :: n_mu, i_r, l, k
        REAL*8, ALLOCATABLE :: mu_q(:), w_q(:), p(:)
        REAL*8 :: mu, rho_val, eta, factor
        REAL*8, PARAMETER :: pi_proj = 3.14159265358979323846D0

        n_mu = MAX(4 * (BASIS_LMAX + 1), 40)
        ALLOCATE(mu_q(n_mu), w_q(n_mu), p(0:BASIS_LMAX))
        CALL gauss_legendre_nodes_weights(n_mu, mu_q, w_q)

        ! eta = 1 - 1/q^2  so the density reads exp(-r/s0 * sqrt(1 - eta*mu^2))
        eta = 1.0D0 - 1.0D0 / (q * q)

        BASIS_RHO_L_GRID = 0.0D0
        DO i_r = 1, BASIS_NR
            DO k = 1, n_mu
                mu = mu_q(k)
                rho_val = rho0 * EXP(-BASIS_R_GRID(i_r) / s0 * SQRT(MAX(1.0D0 - eta*mu*mu, 0.0D0)))
                CALL legendre_p_all_axisymmetric(BASIS_LMAX, mu, p)
                DO l = 0, BASIS_LMAX, 2
                    factor = (2*l + 1) * 0.5D0 * w_q(k)
                    BASIS_RHO_L_GRID(l, i_r) = BASIS_RHO_L_GRID(l, i_r) + factor * rho_val * p(l)
                END DO
            END DO
        END DO

        DEALLOCATE(mu_q, w_q, p)
    END SUBROUTINE project_exponential_oblate_halo
    
    SUBROUTINE exponential_oblate_halo(params, N, x, y, z, ax, ay, az, phi)
        ! Axisymmetric exponential oblate halo:
        !   rho(R,z) = rho0 * exp(-1/s0 * sqrt(R^2 + z^2/q^2))
        ! params = [rho0, s0, q]
        !
        ! On first call (or after clear): if no grid has been set up, auto-
        ! initializes with defaults (lmax=20, 100 log-spaced points 1e-4..1e3).
        ! Then projects rho onto Legendre modes and computes potential tables.
        ! Subsequent calls skip projection and go straight to evaluation.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN),  DIMENSION(4) :: params
        REAL*8, INTENT(IN),  DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi
        REAL*8 :: rho0, s0, q

        rho0 = params(1)
        s0   = params(2)
        q    = params(3)

        IF (.NOT. BASIS_EXPANSION_INITIALIZED) THEN
            IF (.NOT. BASIS_GRID_SET) THEN
                CALL default_init_basis_expansion()
            END IF
            CALL project_exponential_oblate_halo(rho0, s0, q)
            CALL compute_phi_tables_from_rho()
            BASIS_EXPANSION_INITIALIZED = .TRUE.
        END IF

        CALL axisymmetricbasisexpansion_eval(N, x, y, z, ax, ay, az, phi)
    END SUBROUTINE exponential_oblate_halo

    SUBROUTINE project_ibata2024halo(rho0, r0, rt, q, gamma, beta)
        ! Project
        ! rho(s) = rho0 * (s/r0)^(-gamma) * (1 + s/r0)^(gamma-beta) * exp(-(s/rt)^2)
        ! with s = r * sqrt(1 - eta*mu^2), eta = 1 - 1/q^2
        ! onto even-l Legendre modes, filling BASIS_RHO_L_GRID.
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: rho0, r0, rt, q, gamma, beta
        INTEGER :: n_mu, i_r, l, k
        REAL*8, ALLOCATABLE :: mu_q(:), w_q(:), p(:)
        REAL*8 :: mu, rho_val, eta, factor, s, x
        REAL*8, PARAMETER :: s_floor = 1.0D-12

        n_mu = MAX(4 * (BASIS_LMAX + 1), 40)
        ALLOCATE(mu_q(n_mu), w_q(n_mu), p(0:BASIS_LMAX))
        CALL gauss_legendre_nodes_weights(n_mu, mu_q, w_q)

        eta = 1.0D0 - 1.0D0 / (q * q)

        BASIS_RHO_L_GRID = 0.0D0
        DO i_r = 1, BASIS_NR
            DO k = 1, n_mu
                mu = mu_q(k)

                s = BASIS_R_GRID(i_r) * SQRT(MAX(1.0D0 - eta*mu*mu, 0.0D0))
                x = MAX(s / r0, s_floor)

                rho_val = rho0 * x**(-gamma) * (1.0D0 + x)**(gamma - beta) * EXP(-(s/rt)**2)

                CALL legendre_p_all_axisymmetric(BASIS_LMAX, mu, p)
                DO l = 0, BASIS_LMAX, 2
                    factor = (2*l + 1) * 0.5D0 * w_q(k)
                    BASIS_RHO_L_GRID(l, i_r) = BASIS_RHO_L_GRID(l, i_r) + factor * rho_val * p(l)
                END DO
            END DO
        END DO

        DEALLOCATE(mu_q, w_q, p)
    END SUBROUTINE project_ibata2024halo    

    SUBROUTINE ibata2024halo(params, N, x, y, z, ax, ay, az, phi)
        ! Axisymmetric Ibata-like halo:
        ! rho(s) = rho0 * (s/r0)^(-gamma) * (1 + s/r0)^(gamma-beta) * exp(-(s/rt)^2)
        ! s = sqrt(R^2 + z^2/q^2), R^2 = x^2 + y^2
        !
        ! params = [rho0, r0, rt, q, gamma, beta]
        !
        ! Same initialization logic as exponential_oblate_halo.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN),  DIMENSION(7) :: params
        REAL*8, INTENT(IN),  DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi
        REAL*8 :: rho0, r0, rt, q, gamma, beta

        rho0  = params(1)
        r0    = params(2)
        rt    = params(3)
        q     = params(4)
        gamma = params(5)
        beta  = params(6)

        IF (.NOT. BASIS_EXPANSION_INITIALIZED) THEN
            IF (.NOT. BASIS_GRID_SET) THEN
                CALL default_init_basis_expansion()
            END IF
            CALL project_ibata2024halo(rho0, r0, rt, q, gamma, beta)
            CALL compute_phi_tables_from_rho()
            BASIS_EXPANSION_INITIALIZED = .TRUE.
        END IF

        CALL axisymmetricbasisexpansion_eval(N, x, y, z, ax, ay, az, phi)
    END SUBROUTINE ibata2024halo

    ! =======================================================================
    ! =======================================================================
    ! =======================================================================
    ! BFE (COMPOSITE BASIS INFRASTRUCTURE)
    ! Multi-component basis expansion setup, component addition, and evaluation
    ! =======================================================================
    ! =======================================================================
    ! =======================================================================

    SUBROUTINE clearaxisymmetriccompositebasisexpansion()
        IMPLICIT NONE

        IF (ALLOCATED(COMPOSITE_KIND))            DEALLOCATE(COMPOSITE_KIND)
        IF (ALLOCATED(COMPOSITE_READY))           DEALLOCATE(COMPOSITE_READY)
        IF (ALLOCATED(COMPOSITE_R_GRID))          DEALLOCATE(COMPOSITE_R_GRID)
        IF (ALLOCATED(COMPOSITE_RHO_L_GRID))      DEALLOCATE(COMPOSITE_RHO_L_GRID)
        IF (ALLOCATED(COMPOSITE_PHI_L_GRID))      DEALLOCATE(COMPOSITE_PHI_L_GRID)
        IF (ALLOCATED(COMPOSITE_DPHI_L_DR_GRID))  DEALLOCATE(COMPOSITE_DPHI_L_DR_GRID)
        IF (ALLOCATED(COMPOSITE_BESSEL_NPARAMS))    DEALLOCATE(COMPOSITE_BESSEL_NPARAMS)
        IF (ALLOCATED(COMPOSITE_BESSEL_PARAMS))      DEALLOCATE(COMPOSITE_BESSEL_PARAMS)
        IF (ALLOCATED(COMPOSITE_DISK_TABLE_META))    DEALLOCATE(COMPOSITE_DISK_TABLE_META)
        IF (ALLOCATED(COMPOSITE_DISK_TABLE_PHI))     DEALLOCATE(COMPOSITE_DISK_TABLE_PHI)
        IF (ALLOCATED(COMPOSITE_DISK_TABLE_DPHI_DR)) DEALLOCATE(COMPOSITE_DISK_TABLE_DPHI_DR)
        IF (ALLOCATED(COMPOSITE_DISK_TABLE_DPHI_DZ)) DEALLOCATE(COMPOSITE_DISK_TABLE_DPHI_DZ)
        IF (ALLOCATED(COMPOSITE_DISK_TABLE_D2PHI_DRDZ)) DEALLOCATE(COMPOSITE_DISK_TABLE_D2PHI_DRDZ)

        COMPOSITE_BASIS_GRID_SET = .FALSE.
        COMPOSITE_BASIS_FINALIZED = .FALSE.
        COMPOSITE_NCOMP = 0
        COMPOSITE_LMAX = -1
        COMPOSITE_NR = -1
    END SUBROUTINE clearaxisymmetriccompositebasisexpansion

    SUBROUTINE initaxisymmetriccompositebasisexpansion(lmax, nr, r_grid, ncomp)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: lmax, nr, ncomp
        REAL*8, INTENT(IN), DIMENSION(nr) :: r_grid

        CALL clearaxisymmetriccompositebasisexpansion()

        ALLOCATE(COMPOSITE_KIND(ncomp), COMPOSITE_READY(ncomp))
        ALLOCATE(COMPOSITE_BESSEL_NPARAMS(ncomp))
        ALLOCATE(COMPOSITE_R_GRID(nr))
        ALLOCATE(COMPOSITE_RHO_L_GRID(0:lmax, nr, ncomp))
        ALLOCATE(COMPOSITE_PHI_L_GRID(0:lmax, nr, ncomp))
        ALLOCATE(COMPOSITE_DPHI_L_DR_GRID(0:lmax, nr, ncomp))
        ALLOCATE(COMPOSITE_BESSEL_PARAMS(COMPOSITE_BESSEL_MAX_PARAMS, ncomp))
        ALLOCATE(COMPOSITE_DISK_TABLE_META(3, ncomp))
        ALLOCATE(COMPOSITE_DISK_TABLE_PHI(COMPOSITE_DISK_TABLE_NR, COMPOSITE_DISK_TABLE_NZ, ncomp))
        ALLOCATE(COMPOSITE_DISK_TABLE_DPHI_DR(COMPOSITE_DISK_TABLE_NR, COMPOSITE_DISK_TABLE_NZ, ncomp))
        ALLOCATE(COMPOSITE_DISK_TABLE_DPHI_DZ(COMPOSITE_DISK_TABLE_NR, COMPOSITE_DISK_TABLE_NZ, ncomp))
        ALLOCATE(COMPOSITE_DISK_TABLE_D2PHI_DRDZ(COMPOSITE_DISK_TABLE_NR, COMPOSITE_DISK_TABLE_NZ, ncomp))

        COMPOSITE_KIND = BFE_KIND_NONE
        COMPOSITE_READY = .FALSE.
        COMPOSITE_BESSEL_NPARAMS = 0
        COMPOSITE_RHO_L_GRID = 0.0D0
        COMPOSITE_PHI_L_GRID = 0.0D0
        COMPOSITE_DPHI_L_DR_GRID = 0.0D0
        COMPOSITE_BESSEL_PARAMS  = 0.0D0
        COMPOSITE_DISK_TABLE_META = 0.0D0
        COMPOSITE_DISK_TABLE_PHI  = 0.0D0
        COMPOSITE_DISK_TABLE_DPHI_DR = 0.0D0
        COMPOSITE_DISK_TABLE_DPHI_DZ = 0.0D0
        COMPOSITE_DISK_TABLE_D2PHI_DRDZ = 0.0D0

        COMPOSITE_LMAX = lmax
        COMPOSITE_NR = nr
        COMPOSITE_NCOMP = ncomp
        COMPOSITE_R_GRID = r_grid
        COMPOSITE_BASIS_GRID_SET = .TRUE.
        COMPOSITE_BASIS_FINALIZED = .FALSE.
    END SUBROUTINE initaxisymmetriccompositebasisexpansion

    SUBROUTINE addcompositeexponentialoblate(component_index, rho0, s0, q)
        ! Legendre basis path: preferred for spherical or near-spherical components.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: component_index
        REAL*8, INTENT(IN) :: rho0, s0, q
        INTEGER :: n_mu, i_r, l, k
        REAL*8, ALLOCATABLE :: mu_q(:), w_q(:), p(:)
        REAL*8 :: mu, rho_val, eta, factor

        IF (.NOT. COMPOSITE_BASIS_GRID_SET) STOP "initaxisymmetriccompositebasisexpansion must be called before addcompositeexponentialoblate"
        IF (component_index < 1 .OR. component_index > COMPOSITE_NCOMP) STOP "invalid composite component_index"

        n_mu = MAX(4 * (COMPOSITE_LMAX + 1), 40)
        ALLOCATE(mu_q(n_mu), w_q(n_mu), p(0:COMPOSITE_LMAX))
        CALL gauss_legendre_nodes_weights(n_mu, mu_q, w_q)

        eta = 1.0D0 - 1.0D0 / (q * q)
        COMPOSITE_RHO_L_GRID(:,:,component_index) = 0.0D0
        DO i_r = 1, COMPOSITE_NR
            DO k = 1, n_mu
                mu = mu_q(k)
                rho_val = rho0 * EXP(-COMPOSITE_R_GRID(i_r) / s0 * SQRT(MAX(1.0D0 - eta*mu*mu, 0.0D0)))
                CALL legendre_p_all_axisymmetric(COMPOSITE_LMAX, mu, p)
                DO l = 0, COMPOSITE_LMAX, 2
                    factor = (2*l + 1) * 0.5D0 * w_q(k)
                    COMPOSITE_RHO_L_GRID(l, i_r, component_index) = COMPOSITE_RHO_L_GRID(l, i_r, component_index) + factor * rho_val * p(l)
                END DO
            END DO
        END DO
        DEALLOCATE(mu_q, w_q, p)

        CALL compute_phi_tables_from_rho_component(COMPOSITE_R_GRID, &
            COMPOSITE_RHO_L_GRID(:,:,component_index), COMPOSITE_PHI_L_GRID(:,:,component_index), &
            COMPOSITE_DPHI_L_DR_GRID(:,:,component_index))

        COMPOSITE_KIND(component_index) = BFE_KIND_LEGENDRE
        COMPOSITE_READY(component_index) = .TRUE.
        COMPOSITE_BASIS_FINALIZED = .FALSE.
    END SUBROUTINE addcompositeexponentialoblate

    SUBROUTINE addcompositeibata2024halo(component_index, rho0, r0, rt, q, gamma, beta)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: component_index
        REAL*8, INTENT(IN) :: rho0, r0, rt, q, gamma, beta
        INTEGER :: n_mu, i_r, l, k
        REAL*8, ALLOCATABLE :: mu_q(:), w_q(:), p(:)
        REAL*8 :: mu, rho_val, eta, factor, s, x
        REAL*8, PARAMETER :: s_floor = 1.0D-12

        IF (.NOT. COMPOSITE_BASIS_GRID_SET) STOP "initaxisymmetriccompositebasisexpansion must be called before addcompositeibata2024"
        IF (component_index < 1 .OR. component_index > COMPOSITE_NCOMP) STOP "invalid composite component_index"

        n_mu = MAX(4 * (COMPOSITE_LMAX + 1), 40)
        ALLOCATE(mu_q(n_mu), w_q(n_mu), p(0:COMPOSITE_LMAX))
        CALL gauss_legendre_nodes_weights(n_mu, mu_q, w_q)

        eta = 1.0D0 - 1.0D0 / (q * q)
        COMPOSITE_RHO_L_GRID(:,:,component_index) = 0.0D0
        DO i_r = 1, COMPOSITE_NR
            DO k = 1, n_mu
                mu = mu_q(k)
                s = COMPOSITE_R_GRID(i_r) * SQRT(MAX(1.0D0 - eta*mu*mu, 0.0D0))
                x = MAX(s / r0, s_floor)
                rho_val = rho0 * x**(-gamma) * (1.0D0 + x)**(gamma - beta) * EXP(-(s/rt)**2)
                CALL legendre_p_all_axisymmetric(COMPOSITE_LMAX, mu, p)
                DO l = 0, COMPOSITE_LMAX, 2
                    factor = (2*l + 1) * 0.5D0 * w_q(k)
                    COMPOSITE_RHO_L_GRID(l, i_r, component_index) = COMPOSITE_RHO_L_GRID(l, i_r, component_index) + factor * rho_val * p(l)
                END DO
            END DO
        END DO
        DEALLOCATE(mu_q, w_q, p)

        CALL compute_phi_tables_from_rho_component(COMPOSITE_R_GRID, &
            COMPOSITE_RHO_L_GRID(:,:,component_index), COMPOSITE_PHI_L_GRID(:,:,component_index), &
            COMPOSITE_DPHI_L_DR_GRID(:,:,component_index))

        COMPOSITE_KIND(component_index) = BFE_KIND_LEGENDRE
        COMPOSITE_READY(component_index) = .TRUE.
        COMPOSITE_BASIS_FINALIZED = .FALSE.
    END SUBROUTINE addcompositeibata2024halo

    SUBROUTINE addcompositebesselcomponent(component_index, params, nparams)
        ! Generic Bessel component registration.
        ! params is evaluator-specific and interpreted by the selected
        ! Bessel application evaluator.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: component_index
        INTEGER, INTENT(IN) :: nparams
        REAL*8, INTENT(IN), DIMENSION(nparams) :: params

        IF (.NOT. COMPOSITE_BASIS_GRID_SET) STOP "initaxisymmetriccompositebasisexpansion must be called before addcompositebesselcomponent"
        IF (component_index < 1 .OR. component_index > COMPOSITE_NCOMP) STOP "invalid composite component_index"
        IF (nparams < 1) STOP "nparams must be >= 1 in addcompositebesselcomponent"
        IF (nparams > COMPOSITE_BESSEL_MAX_PARAMS) STOP "nparams exceeds COMPOSITE_BESSEL_MAX_PARAMS"

        COMPOSITE_BESSEL_PARAMS(:, component_index) = 0.0D0
        COMPOSITE_BESSEL_PARAMS(1:nparams, component_index) = params(1:nparams)
        COMPOSITE_BESSEL_NPARAMS(component_index) = nparams
        COMPOSITE_KIND(component_index) = BFE_KIND_BESSEL_DISK
        COMPOSITE_READY(component_index) = .TRUE.
        COMPOSITE_BASIS_FINALIZED = .FALSE.
    END SUBROUTINE addcompositebesselcomponent

    SUBROUTINE addcompositebesselexponentialdisk(component_index, sigma0, hR, hZ)
        ! Register and precompute an exponential-disk component using the
        ! tabulated cylindrical force representation (BFE_KIND_DISK_TABLE).
        ! The one-time Bessel/Hankel quadrature is performed here at setup time;
        ! runtime evaluation uses bilinear interpolation on the stored table.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: component_index
        REAL*8, INTENT(IN) :: sigma0, hR, hZ
        REAL*8, DIMENSION(4) :: params_disk

        IF (.NOT. COMPOSITE_BASIS_GRID_SET) STOP "initaxisymmetriccompositebasisexpansion must be called before addcompositebesselexponentialdisk"
        IF (component_index < 1 .OR. component_index > COMPOSITE_NCOMP) STOP "invalid composite component_index"
        IF (hR <= 0.0D0 .OR. hZ <= 0.0D0) STOP "hR and hZ must be positive"

        ! Retain raw params for the reference direct-quadrature path.
        params_disk(1) = sigma0
        params_disk(2) = hR
        params_disk(3) = hZ
        COMPOSITE_BESSEL_PARAMS(:, component_index) = 0.0D0
        COMPOSITE_BESSEL_PARAMS(1:4, component_index) = params_disk
        COMPOSITE_BESSEL_NPARAMS(component_index) = 4

        ! Build the 2D cylindrical table (expensive offline step).
        CALL build_exponential_disk_table(component_index, sigma0, hR, hZ)

        COMPOSITE_KIND(component_index) = BFE_KIND_DISK_TABLE
        COMPOSITE_READY(component_index) = .TRUE.
        COMPOSITE_BASIS_FINALIZED = .FALSE.
    END SUBROUTINE addcompositebesselexponentialdisk

    SUBROUTINE finalizeaxisymmetriccompositebasisexpansion()
        IMPLICIT NONE
        INTEGER :: i

        IF (.NOT. COMPOSITE_BASIS_GRID_SET) STOP "initaxisymmetriccompositebasisexpansion must be called before finalizeaxisymmetriccompositebasisexpansion"
        DO i = 1, COMPOSITE_NCOMP
            IF (.NOT. COMPOSITE_READY(i)) STOP "all composite components must be configured before finalizeaxisymmetriccompositebasisexpansion"
        END DO
        COMPOSITE_BASIS_FINALIZED = .TRUE.
    END SUBROUTINE finalizeaxisymmetriccompositebasisexpansion

    SUBROUTINE axisymmetriccompositebasispotential_dispatch(params, N, x, y, z, ax, ay, az, phi)
        ! Compatibility bridge for simulator's generic params-based pointer API.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(*) :: params
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi

        ! Keep params consumed so compilers do not warn in this wrapper.
        IF (params(1) /= params(1)) THEN
            ax = 0.0D0
            ay = 0.0D0
            az = 0.0D0
            phi = 0.0D0
            RETURN
        END IF

        CALL axisymmetriccompositebasispotential(N, x, y, z, ax, ay, az, phi)
    END SUBROUTINE axisymmetriccompositebasispotential_dispatch

    SUBROUTINE axisymmetriccompositebasispotential(N, x, y, z, ax, ay, az, phi)
        ! Evaluate all configured composite basis components and sum their forces.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi
        REAL*8, ALLOCATABLE :: ax_comp(:), ay_comp(:), az_comp(:), phi_comp(:)
        INTEGER :: i, nparams_bessel

        ax = 0.0D0
        ay = 0.0D0
        az = 0.0D0
        phi = 0.0D0

        IF (.NOT. COMPOSITE_BASIS_FINALIZED) STOP "finalizeaxisymmetriccompositebasisexpansion must be called before axisymmetriccompositebasispotential"

        ALLOCATE(ax_comp(N), ay_comp(N), az_comp(N), phi_comp(N))
        DO i = 1, COMPOSITE_NCOMP
            ax_comp = 0.0D0
            ay_comp = 0.0D0
            az_comp = 0.0D0
            phi_comp = 0.0D0

            SELECT CASE (COMPOSITE_KIND(i))
            CASE (BFE_KIND_LEGENDRE)
                CALL axisymmetricbasisexpansion_eval_component(N, x, y, z, COMPOSITE_R_GRID, &
                    COMPOSITE_PHI_L_GRID(:,:,i), COMPOSITE_DPHI_L_DR_GRID(:,:,i), &
                    ax_comp, ay_comp, az_comp, phi_comp)
            CASE (BFE_KIND_BESSEL_DISK)
                nparams_bessel = COMPOSITE_BESSEL_NPARAMS(i)
                IF (nparams_bessel < 1) STOP "invalid bessel parameter count in axisymmetriccompositebasispotential"
                CALL bessel_eval_component(exponential_disk_bessel_eval_component, &
                    COMPOSITE_BESSEL_PARAMS(1:nparams_bessel, i), N, &
                    x, y, z, ax_comp, ay_comp, az_comp, phi_comp)
            CASE (BFE_KIND_DISK_TABLE)
                CALL disk_table_eval_component(i, N, x, y, z, ax_comp, ay_comp, az_comp, phi_comp)
            CASE DEFAULT
                STOP "unknown component kind in axisymmetriccompositebasispotential"
            END SELECT

            ax = ax + ax_comp
            ay = ay + ay_comp
            az = az + az_comp
            phi = phi + phi_comp
        END DO
        DEALLOCATE(ax_comp, ay_comp, az_comp, phi_comp)
    END SUBROUTINE axisymmetriccompositebasispotential

    SUBROUTINE bessel_eval_component(component_evaluator, params, N, x, y, z, ax, ay, az, phi)
        ! Generic Bessel-component dispatcher.
        ! Method: Bessel/Hankel representation used to solve Poisson for
        ! flattened axisymmetric components.
        ! Application: the specific density profile is implemented by the
        ! passed component_evaluator.
        IMPLICIT NONE
        PROCEDURE(bessel_component_evaluator_interface) :: component_evaluator
        REAL*8, INTENT(IN), DIMENSION(*) :: params
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi

        CALL component_evaluator(params, N, x, y, z, ax, ay, az, phi)
    END SUBROUTINE bessel_eval_component

end module gravity




