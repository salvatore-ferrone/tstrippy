MODULE gravity
    USE sphericalharmonicsbfe, ONLY: BASIS_GRID_SET, BASIS_EXPANSION_INITIALIZED, &
                                     sh_set_basis_g => setsphericalharmonicbasisgravityconstant, &
                                     sh_clear_basis => clearsphericalharmonicbasis, &
                                     sh_init_basis => initsphericalharmonicbasis, &
                                     sh_default_init_basis => defaultinitsphericalharmonicbasis, &
                                     sh_init_component_phi => initsphericalharmoniccomponentphi, &
                                     sh_store_component_phi => storesphericalharmoniccomponentphi, &
                                     sh_load_component_phi => loadsphericalharmoniccomponentphi, &
                                     sh_project_density => project_axisym_density_generic, &
                                     sh_compute_phi_tables => compute_phi_tables_from_rho, &
                                     sh_eval_force => sphericalharmonicbasisforce, &
                                     sh_eval_potential => sphericalharmonicbasispotential
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

    PUBLIC :: ibata2024halo_density
    PRIVATE :: count_sh_components, sh_slot_for_component

CONTAINS

    SUBROUTINE cleargravity()
        IMPLICIT NONE
        GRAVITY_G = GRAVITY_G_DEFAULT
        GRAVITY_G_IS_DEFAULT = .TRUE.
        GRAVITY_FINALIZED = .FALSE.
        GRAVITY_NCOMP = 0
        GRAVITY_KIND = GRAVITY_KIND_NONE
        GRAVITY_PARAMS = 0.0D0
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
        GRAVITY_PARAMS(1:nparams, GRAVITY_NCOMP) = params(1:nparams)
    END SUBROUTINE addgravitycomponent

    SUBROUTINE finalizegravity()
        IMPLICIT NONE
        INTEGER :: i, n_sh, i_sh, i_sh_slot
        IF (GRAVITY_NCOMP < 1) THEN
            WRITE(*,'(A)') "WARNING: finalizegravity: no components registered"
            RETURN
        END IF

        n_sh = 0
        i_sh = 0
        DO i = 1, GRAVITY_NCOMP
            IF (GRAVITY_KIND(i) == GRAVITY_KIND_NONE) THEN
                WRITE(*,'(A)') "WARNING: finalizegravity: component slot is uninitialized"
                RETURN
            END IF
            IF (GRAVITY_KIND(i) == GRAVITY_KIND_EXPONENTIALOBLATEHALO .OR. &
                GRAVITY_KIND(i) == GRAVITY_KIND_IBATA2024HALO) THEN
                n_sh = n_sh + 1
                i_sh = i
            END IF
        END DO

        ! Microstep: eager SH table build in finalize for the single-SH-component case.
        IF (n_sh == 1) THEN
            IF (.NOT. BASIS_GRID_SET) CALL sh_default_init_basis()
            SELECT CASE (GRAVITY_KIND(i_sh))
            CASE (GRAVITY_KIND_EXPONENTIALOBLATEHALO)
                CALL sh_project_density(GRAVITY_PARAMS(1:3, i_sh), exponentialoblatehalo_density)
            CASE (GRAVITY_KIND_IBATA2024HALO)
                CALL sh_project_density(GRAVITY_PARAMS(1:6, i_sh), ibata2024halo_density)
            END SELECT
            CALL sh_compute_phi_tables()
        ELSE IF (n_sh > 1) THEN
            IF (.NOT. BASIS_GRID_SET) CALL sh_default_init_basis()
            CALL sh_init_component_phi(n_sh)

            i_sh_slot = 0
            DO i = 1, GRAVITY_NCOMP
                SELECT CASE (GRAVITY_KIND(i))
                CASE (GRAVITY_KIND_EXPONENTIALOBLATEHALO)
                    i_sh_slot = i_sh_slot + 1
                    CALL sh_project_density(GRAVITY_PARAMS(1:3, i), exponentialoblatehalo_density)
                    CALL sh_compute_phi_tables()
                    CALL sh_store_component_phi(i_sh_slot)
                CASE (GRAVITY_KIND_IBATA2024HALO)
                    i_sh_slot = i_sh_slot + 1
                    CALL sh_project_density(GRAVITY_PARAMS(1:6, i), ibata2024halo_density)
                    CALL sh_compute_phi_tables()
                    CALL sh_store_component_phi(i_sh_slot)
                END SELECT
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
        REAL*8, DIMENSION(n) :: ax_c, ay_c, az_c
        INTEGER :: i, n_sh

        ax_comp = 0.0D0
        ay_comp = 0.0D0
        az_comp = 0.0D0

        IF (.NOT. GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: evaluategravityforcecomponents: call finalizegravity first"
            RETURN
        END IF

        n_sh = count_sh_components()

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
            CASE (GRAVITY_KIND_EXPONENTIALOBLATEHALO)
                IF (n_sh > 1) THEN
                    CALL sh_load_component_phi(sh_slot_for_component(i))
                ELSE IF (.NOT. BASIS_EXPANSION_INITIALIZED) THEN
                    IF (.NOT. BASIS_GRID_SET) CALL sh_default_init_basis()
                    CALL sh_project_density(GRAVITY_PARAMS(1:3, i), exponentialoblatehalo_density)
                    CALL sh_compute_phi_tables()
                END IF
                CALL sh_eval_force(n, x, y, z, ax_c, ay_c, az_c)
                ax_comp(i,:) = ax_c
                ay_comp(i,:) = ay_c
                az_comp(i,:) = az_c
            CASE (GRAVITY_KIND_IBATA2024HALO)
                IF (n_sh > 1) THEN
                    CALL sh_load_component_phi(sh_slot_for_component(i))
                ELSE IF (.NOT. BASIS_EXPANSION_INITIALIZED) THEN
                    IF (.NOT. BASIS_GRID_SET) CALL sh_default_init_basis()
                    CALL sh_project_density(GRAVITY_PARAMS(1:6, i), ibata2024halo_density)
                    CALL sh_compute_phi_tables()
                END IF
                CALL sh_eval_force(n, x, y, z, ax_c, ay_c, az_c)
                ax_comp(i,:) = ax_c
                ay_comp(i,:) = ay_c
                az_comp(i,:) = az_c
            END SELECT
        END DO
    END SUBROUTINE evaluategravityforcecomponents

    INTEGER FUNCTION count_sh_components()
        IMPLICIT NONE
        INTEGER :: i

        count_sh_components = 0
        DO i = 1, GRAVITY_NCOMP
            IF (GRAVITY_KIND(i) == GRAVITY_KIND_EXPONENTIALOBLATEHALO .OR. &
                GRAVITY_KIND(i) == GRAVITY_KIND_IBATA2024HALO) THEN
                count_sh_components = count_sh_components + 1
            END IF
        END DO
    END FUNCTION count_sh_components

    INTEGER FUNCTION sh_slot_for_component(i_comp)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: i_comp
        INTEGER :: i

        sh_slot_for_component = 0
        IF (i_comp < 1 .OR. i_comp > GRAVITY_NCOMP) RETURN

        DO i = 1, i_comp
            IF (GRAVITY_KIND(i) == GRAVITY_KIND_EXPONENTIALOBLATEHALO .OR. &
                GRAVITY_KIND(i) == GRAVITY_KIND_IBATA2024HALO) THEN
                sh_slot_for_component = sh_slot_for_component + 1
            END IF
        END DO
    END FUNCTION sh_slot_for_component

    SUBROUTINE evaluategravityforces(n, x, y, z, ax, ay, az)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: ax, ay, az
        REAL*8, DIMENSION(n,3) :: force_tmp
        REAL*8, DIMENSION(n) :: ax_c, ay_c, az_c
        INTEGER :: i, n_sh

        ax = 0.0D0
        ay = 0.0D0
        az = 0.0D0

        IF (.NOT. GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: evaluategravityforces: call finalizegravity first"
            RETURN
        END IF

        n_sh = count_sh_components()

        DO i = 1, GRAVITY_NCOMP
            SELECT CASE (GRAVITY_KIND(i))
            CASE (GRAVITY_KIND_PLUMMER)
                CALL plummer_force(GRAVITY_PARAMS(1:2, i), n, x, y, z, force_tmp)
                ax = ax + force_tmp(:,1)
                ay = ay + force_tmp(:,2)
                az = az + force_tmp(:,3)
            CASE (GRAVITY_KIND_HERNQUIST)
                CALL hernquist_force(GRAVITY_PARAMS(1:2, i), n, x, y, z, force_tmp)
                ax = ax + force_tmp(:,1)
                ay = ay + force_tmp(:,2)
                az = az + force_tmp(:,3)
            CASE (GRAVITY_KIND_EXPONENTIALOBLATEHALO)
                IF (n_sh > 1) THEN
                    CALL sh_load_component_phi(sh_slot_for_component(i))
                ELSE IF (.NOT. BASIS_EXPANSION_INITIALIZED) THEN
                    IF (.NOT. BASIS_GRID_SET) CALL sh_default_init_basis()
                    CALL sh_project_density(GRAVITY_PARAMS(1:3, i), exponentialoblatehalo_density)
                    CALL sh_compute_phi_tables()
                END IF
                CALL sh_eval_force(n, x, y, z, ax_c, ay_c, az_c)
                ax = ax + ax_c
                ay = ay + ay_c
                az = az + az_c
            CASE (GRAVITY_KIND_IBATA2024HALO)
                IF (n_sh > 1) THEN
                    CALL sh_load_component_phi(sh_slot_for_component(i))
                ELSE IF (.NOT. BASIS_EXPANSION_INITIALIZED) THEN
                    IF (.NOT. BASIS_GRID_SET) CALL sh_default_init_basis()
                    CALL sh_project_density(GRAVITY_PARAMS(1:6, i), ibata2024halo_density)
                    CALL sh_compute_phi_tables()
                END IF
                CALL sh_eval_force(n, x, y, z, ax_c, ay_c, az_c)
                ax = ax + ax_c
                ay = ay + ay_c
                az = az + az_c
            END SELECT
        END DO
    END SUBROUTINE evaluategravityforces

    SUBROUTINE evaluategravitypotential(n, x, y, z, phi)
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: phi
        REAL*8, DIMENSION(n) :: phi_c
        INTEGER :: i, n_sh

        phi = 0.0D0

        IF (.NOT. GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: evaluategravitypotential: call finalizegravity first"
            RETURN
        END IF

        n_sh = count_sh_components()

        DO i = 1, GRAVITY_NCOMP
            SELECT CASE (GRAVITY_KIND(i))
            CASE (GRAVITY_KIND_PLUMMER)
                CALL plummer_potential(GRAVITY_PARAMS(1:2, i), n, x, y, z, phi_c)
                phi = phi + phi_c
            CASE (GRAVITY_KIND_HERNQUIST)
                CALL hernquist_potential(GRAVITY_PARAMS(1:2, i), n, x, y, z, phi_c)
                phi = phi + phi_c
            CASE (GRAVITY_KIND_EXPONENTIALOBLATEHALO)
                IF (n_sh > 1) THEN
                    CALL sh_load_component_phi(sh_slot_for_component(i))
                ELSE IF (.NOT. BASIS_EXPANSION_INITIALIZED) THEN
                    IF (.NOT. BASIS_GRID_SET) CALL sh_default_init_basis()
                    CALL sh_project_density(GRAVITY_PARAMS(1:3, i), exponentialoblatehalo_density)
                    CALL sh_compute_phi_tables()
                END IF
                CALL sh_eval_potential(n, x, y, z, phi_c)
                phi = phi + phi_c
            CASE (GRAVITY_KIND_IBATA2024HALO)
                IF (n_sh > 1) THEN
                    CALL sh_load_component_phi(sh_slot_for_component(i))
                ELSE IF (.NOT. BASIS_EXPANSION_INITIALIZED) THEN
                    IF (.NOT. BASIS_GRID_SET) CALL sh_default_init_basis()
                    CALL sh_project_density(GRAVITY_PARAMS(1:6, i), ibata2024halo_density)
                    CALL sh_compute_phi_tables()
                END IF
                CALL sh_eval_potential(n, x, y, z, phi_c)
                phi = phi + phi_c
            END SELECT
        END DO
    END SUBROUTINE evaluategravitypotential

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! ANALYTICAL POTENTIAL MODELS !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! DENSITY ONLY PROFILES !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

    SUBROUTINE ibata2024halo_density(params, n, x, y, z, rho)
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
    END SUBROUTINE ibata2024halo_density    

END MODULE gravity
