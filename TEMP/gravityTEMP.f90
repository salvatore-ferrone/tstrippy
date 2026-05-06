MODULE gravity
    USE sphericalharmonicsbfe, ONLY: BASIS_GRID_SET, BASIS_EXPANSION_INITIALIZED, BASIS_RHO_L_GRID, &
                                     BASIS_RHO_L_COMPONENT_GRID, BASIS_PHI_L_COMPONENT_GRID, &
                                     BASIS_DPHI_L_DR_COMPONENT_GRID, &
                                     sh_set_basis_g => setsphericalharmonicbasisgravityconstant, &
                                     sh_clear_basis => clearsphericalharmonicbasis, &
                                     sh_init_basis => initsphericalharmonicbasis, &
                                     sh_default_init_basis => defaultinitsphericalharmonicbasis, &
                                     sh_init_component_storage => initsphericalharmoniccomponentstorage, &
                                     sh_project_density_to_table => project_axisymmetric_density_to_table, &
                                     sh_compute_phi_tables => compute_phi_tables_from_rho, &
                                     sh_compute_phi_from_rho_input => compute_phi_tables_from_rho_input, &
                                     sh_eval_force_component => sphericalharmonicbasisforce_component, &
                                     sh_eval_potential_component => sphericalharmonicbasispotential_component
    IMPLICIT NONE

    REAL*8, PARAMETER, PUBLIC :: GRAVITY_G_DEFAULT = 4.30091727D-6
    INTEGER, PARAMETER, PUBLIC :: GRAVITY_MAX_NCOMP = 16
    INTEGER, PARAMETER, PUBLIC :: GRAVITY_MAX_PARAMS = 16

    INTEGER, PARAMETER, PUBLIC :: GRAVITY_KIND_NONE = 0
    INTEGER, PARAMETER, PUBLIC :: GRAVITY_KIND_PLUMMER = 10
    INTEGER, PARAMETER, PUBLIC :: GRAVITY_KIND_HERNQUIST = 11
    INTEGER, PARAMETER, PUBLIC :: GRAVITY_KIND_EXPONENTIALOBLATEHALO = 20
    INTEGER, PARAMETER, PUBLIC :: GRAVITY_KIND_IBATA2024HALO = 21

    REAL*8, PUBLIC :: GRAVITY_G = GRAVITY_G_DEFAULT
    LOGICAL, PUBLIC :: GRAVITY_G_IS_DEFAULT = .TRUE.
    LOGICAL, PUBLIC :: GRAVITY_FINALIZED = .FALSE.
    INTEGER, PUBLIC :: GRAVITY_NCOMP = 0
    INTEGER, DIMENSION(GRAVITY_MAX_NCOMP), PUBLIC :: GRAVITY_KIND = 0
    REAL*8, DIMENSION(GRAVITY_MAX_PARAMS, GRAVITY_MAX_NCOMP), PUBLIC :: GRAVITY_PARAMS = 0.0D0

    INTEGER, DIMENSION(GRAVITY_MAX_NCOMP) :: GRAVITY_SH_SLOT_FOR_COMPONENT = 0

CONTAINS

    SUBROUTINE cleargravity()
        IMPLICIT NONE
        GRAVITY_G = GRAVITY_G_DEFAULT
        GRAVITY_G_IS_DEFAULT = .TRUE.
        GRAVITY_FINALIZED = .FALSE.
        GRAVITY_NCOMP = 0
        GRAVITY_KIND = GRAVITY_KIND_NONE
        GRAVITY_PARAMS = 0.0D0
        GRAVITY_SH_SLOT_FOR_COMPONENT = 0
        CALL sh_clear_basis()
        CALL sh_set_basis_g(GRAVITY_G)
    END SUBROUTINE cleargravity

    SUBROUTINE setgravityconstant(g)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: g
        IF (GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: setgravityconstant: cannot change G after finalizegravity"
            RETURN
        END IF
        IF (g <= 0.0D0) THEN
            WRITE(*,'(A)') "WARNING: setgravityconstant: G must be positive"
            RETURN
        END IF
        GRAVITY_G = g
        GRAVITY_G_IS_DEFAULT = .FALSE.
        CALL sh_set_basis_g(GRAVITY_G)
    END SUBROUTINE setgravityconstant

    SUBROUTINE clearsphericalharmonicbasis()
        IMPLICIT NONE
        CALL sh_clear_basis()
    END SUBROUTINE clearsphericalharmonicbasis

    SUBROUTINE initsphericalharmonicbasis(lmax, nr, r_grid)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: lmax, nr
        REAL*8, INTENT(IN), DIMENSION(nr) :: r_grid
        CALL sh_set_basis_g(GRAVITY_G)
        CALL sh_init_basis(lmax, nr, r_grid)
    END SUBROUTINE initsphericalharmonicbasis

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

        SELECT CASE (TRIM(model_name))
        CASE ("plummer")
            kind_code = GRAVITY_KIND_PLUMMER
            required_params = 2
        CASE ("hernquist")
            kind_code = GRAVITY_KIND_HERNQUIST
            required_params = 2
        CASE ("exponentialoblatehalo")
            kind_code = GRAVITY_KIND_EXPONENTIALOBLATEHALO
            required_params = 3
        CASE ("ibata2024halo")
            kind_code = GRAVITY_KIND_IBATA2024HALO
            required_params = 6
        CASE DEFAULT
            WRITE(*,'(A)') "WARNING: addgravitycomponent: unknown model_name"
            RETURN
        END SELECT

        IF (nparams /= required_params) THEN
            WRITE(*,'(A,I0,A,I0)') "WARNING: addgravitycomponent: expected", required_params, " params, got", nparams
            RETURN
        END IF

        GRAVITY_NCOMP = GRAVITY_NCOMP + 1
        GRAVITY_KIND(GRAVITY_NCOMP) = kind_code
        GRAVITY_PARAMS(:, GRAVITY_NCOMP) = 0.0D0
        GRAVITY_PARAMS(1:nparams, GRAVITY_NCOMP) = params(1:nparams)
    END SUBROUTINE addgravitycomponent

    SUBROUTINE finalizegravity()
        IMPLICIT NONE
        INTEGER :: i, sh_count, sh_slot

        IF (GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: finalizegravity: already finalized"
            RETURN
        END IF
        IF (GRAVITY_NCOMP < 1) THEN
            WRITE(*,'(A)') "WARNING: finalizegravity: no components registered"
            RETURN
        END IF

        GRAVITY_SH_SLOT_FOR_COMPONENT = 0
        sh_count = 0

        DO i = 1, GRAVITY_NCOMP
            IF (GRAVITY_KIND(i) == GRAVITY_KIND_NONE) THEN
                WRITE(*,'(A)') "WARNING: finalizegravity: component slot is uninitialized"
                RETURN
            END IF
            IF (gravity_kind_uses_sh(GRAVITY_KIND(i))) THEN
                sh_count = sh_count + 1
                GRAVITY_SH_SLOT_FOR_COMPONENT(i) = sh_count
            END IF
        END DO

        IF (sh_count > 0) THEN
            IF (.NOT. BASIS_GRID_SET) CALL sh_default_init_basis()
            CALL sh_set_basis_g(GRAVITY_G)
            CALL sh_init_component_storage(sh_count)

            BASIS_RHO_L_GRID = 0.0D0
            DO i = 1, GRAVITY_NCOMP
                sh_slot = GRAVITY_SH_SLOT_FOR_COMPONENT(i)
                IF (sh_slot < 1) CYCLE

                SELECT CASE (GRAVITY_KIND(i))
                CASE (GRAVITY_KIND_EXPONENTIALOBLATEHALO)
                    CALL sh_project_density_to_table(GRAVITY_PARAMS(1:3, i), exponentialoblatehalo_density, &
                                                     BASIS_RHO_L_COMPONENT_GRID(:,:,sh_slot))
                CASE (GRAVITY_KIND_IBATA2024HALO)
                    CALL sh_project_density_to_table(GRAVITY_PARAMS(1:6, i), haloibata2024_density, &
                                                     BASIS_RHO_L_COMPONENT_GRID(:,:,sh_slot))
                CASE DEFAULT
                    CYCLE
                END SELECT

                BASIS_RHO_L_GRID = BASIS_RHO_L_GRID + BASIS_RHO_L_COMPONENT_GRID(:,:,sh_slot)
            END DO

            CALL sh_compute_phi_tables()

            DO sh_slot = 1, sh_count
                CALL sh_compute_phi_from_rho_input(BASIS_RHO_L_COMPONENT_GRID(:,:,sh_slot), &
                                                   BASIS_PHI_L_COMPONENT_GRID(:,:,sh_slot), &
                                                   BASIS_DPHI_L_DR_COMPONENT_GRID(:,:,sh_slot))
            END DO
        END IF

        GRAVITY_FINALIZED = .TRUE.
    END SUBROUTINE finalizegravity

    SUBROUTINE evaluategravityforcecomponents(n, x, y, z, ax_comp, ay_comp, az_comp)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(16, n) :: ax_comp, ay_comp, az_comp
        REAL*8, DIMENSION(n,3) :: force_tmp
        INTEGER :: i, sh_slot

        ax_comp = 0.0D0
        ay_comp = 0.0D0
        az_comp = 0.0D0

        IF (.NOT. GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: evaluategravityforcecomponents: call finalizegravity first"
            RETURN
        END IF

        DO i = 1, GRAVITY_NCOMP
            SELECT CASE (GRAVITY_KIND(i))
            CASE (GRAVITY_KIND_PLUMMER)
                CALL plummer_force(GRAVITY_PARAMS(1:2, i), n, x, y, z, force_tmp)
                ax_comp(i,:) = force_tmp(:,1)
                ay_comp(i,:) = force_tmp(:,2)
                az_comp(i,:) = force_tmp(:,3)
            CASE (GRAVITY_KIND_HERNQUIST)
                CALL hernquist_force(GRAVITY_PARAMS(1:2, i), n, x, y, z, force_tmp)
                ax_comp(i,:) = force_tmp(:,1)
                ay_comp(i,:) = force_tmp(:,2)
                az_comp(i,:) = force_tmp(:,3)
            CASE (GRAVITY_KIND_EXPONENTIALOBLATEHALO, GRAVITY_KIND_IBATA2024HALO)
                sh_slot = GRAVITY_SH_SLOT_FOR_COMPONENT(i)
                IF (sh_slot > 0) THEN
                    CALL sh_eval_force_component(sh_slot, n, x, y, z, ax_comp(i,:), ay_comp(i,:), az_comp(i,:))
                END IF
            END SELECT
        END DO
    END SUBROUTINE evaluategravityforcecomponents

    SUBROUTINE evaluategravityforces(n, x, y, z, ax, ay, az)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: ax, ay, az
        REAL*8 :: ax_comp(16, n), ay_comp(16, n), az_comp(16, n)
        INTEGER :: i

        ax = 0.0D0
        ay = 0.0D0
        az = 0.0D0

        IF (.NOT. GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: evaluategravityforces: call finalizegravity first"
            RETURN
        END IF

        CALL evaluategravityforcecomponents(n, x, y, z, ax_comp, ay_comp, az_comp)

        DO i = 1, GRAVITY_NCOMP
            ax = ax + ax_comp(i,:)
            ay = ay + ay_comp(i,:)
            az = az + az_comp(i,:)
        END DO

    END SUBROUTINE evaluategravityforces

    SUBROUTINE evaluategravitypotential(n, x, y, z, phi)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: phi
        REAL*8, DIMENSION(n) :: phi_c
        INTEGER :: i, sh_slot

        phi = 0.0D0

        IF (.NOT. GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: evaluategravitypotential: call finalizegravity first"
            RETURN
        END IF

        DO i = 1, GRAVITY_NCOMP
            SELECT CASE (GRAVITY_KIND(i))
            CASE (GRAVITY_KIND_PLUMMER)
                CALL plummer_potential(GRAVITY_PARAMS(1:2, i), n, x, y, z, phi_c)
                phi = phi + phi_c
            CASE (GRAVITY_KIND_HERNQUIST)
                CALL hernquist_potential(GRAVITY_PARAMS(1:2, i), n, x, y, z, phi_c)
                phi = phi + phi_c
            CASE (GRAVITY_KIND_EXPONENTIALOBLATEHALO, GRAVITY_KIND_IBATA2024HALO)
                sh_slot = GRAVITY_SH_SLOT_FOR_COMPONENT(i)
                IF (sh_slot > 0) THEN
                    CALL sh_eval_potential_component(sh_slot, n, x, y, z, phi_c)
                    phi = phi + phi_c
                END IF
            END SELECT
        END DO
    END SUBROUTINE evaluategravitypotential

    LOGICAL FUNCTION gravity_kind_uses_sh(kind_code)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: kind_code
        SELECT CASE (kind_code)
        CASE (GRAVITY_KIND_EXPONENTIALOBLATEHALO, GRAVITY_KIND_IBATA2024HALO)
            gravity_kind_uses_sh = .TRUE.
        CASE DEFAULT
            gravity_kind_uses_sh = .FALSE.
        END SELECT
    END FUNCTION gravity_kind_uses_sh

    SUBROUTINE plummer_force(params, n, x, y, z, force)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(2) :: params
        REAL*8, INTENT(OUT), DIMENSION(n,3) :: force
        REAL*8, DIMENSION(n) :: r, amod
        REAL*8 :: m, b

        m = params(1)
        b = params(2)
        r = SQRT(x*x + y*y + z*z)
        amod = -GRAVITY_G*m / (r*r + b*b)**1.5

        force(:,1) = amod*x
        force(:,2) = amod*y
        force(:,3) = amod*z
    END SUBROUTINE plummer_force

    SUBROUTINE plummer_potential(params, n, x, y, z, phi)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(2) :: params
        REAL*8, INTENT(OUT), DIMENSION(n) :: phi
        REAL*8, DIMENSION(n) :: r
        REAL*8 :: m, b

        m = params(1)
        b = params(2)
        r = SQRT(x*x + y*y + z*z)
        phi = -GRAVITY_G*m / SQRT(r*r + b*b)
    END SUBROUTINE plummer_potential

    SUBROUTINE hernquist_force(params, n, x, y, z, force)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(2) :: params
        REAL*8, INTENT(OUT), DIMENSION(n,3) :: force
        REAL*8, DIMENSION(n) :: r, amod
        REAL*8 :: m, a
        REAL*8, PARAMETER :: eps = 1.0D-30

        m = params(1)
        a = params(2)
        r = SQRT(x*x + y*y + z*z)
        amod = -GRAVITY_G*m / (MAX(r, eps) * (r + a)**2)

        force(:,1) = amod*x
        force(:,2) = amod*y
        force(:,3) = amod*z
    END SUBROUTINE hernquist_force

    SUBROUTINE hernquist_potential(params, n, x, y, z, phi)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(2) :: params
        REAL*8, INTENT(OUT), DIMENSION(n) :: phi
        REAL*8, DIMENSION(n) :: r
        REAL*8 :: m, a

        m = params(1)
        a = params(2)
        r = SQRT(x*x + y*y + z*z)
        phi = -GRAVITY_G*m / (r + a)
    END SUBROUTINE hernquist_potential

    SUBROUTINE exponentialoblatehalo_density(params, n, x, y, z, rho)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: rho
        REAL*8 :: rho0, s0, q
        REAL*8, DIMENSION(n) :: r2, s

        rho0 = params(1)
        s0 = params(2)
        q = params(3)

        r2 = x**2 + y**2
        s = SQRT(r2 + (z/q)**2)
        rho = rho0 * EXP(-(1.0D0/s0) * s)
    END SUBROUTINE exponentialoblatehalo_density

    SUBROUTINE haloibata2024_density(params, n, x, y, z, rho)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN),  DIMENSION(:) :: params
        REAL*8, INTENT(IN),  DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: rho
        REAL*8 :: rho0, r0, rt, q, gamma, beta
        REAL*8, DIMENSION(n) :: r2, s, s_safe
        REAL*8, PARAMETER :: eps = 1.0D-30

        rho0  = params(1)
        r0    = params(2)
        rt    = params(3)
        q     = params(4)
        gamma = params(5)
        beta  = params(6)

        r2 = x**2 + y**2
        s = SQRT(r2 + (z/q)**2)
        s_safe = MAX(s, eps)

        rho = rho0 * (s_safe/r0)**(-gamma) * (1.0D0 + s_safe/r0)**(gamma-beta) * EXP(-(s_safe/rt)**2)
    END SUBROUTINE haloibata2024_density

END MODULE gravity
