MODULE gravity
    USE sphericalharmonicsbfe, ONLY: BASIS_GRID_SET, BASIS_EXPANSION_INITIALIZED, &
                                     sh_set_basis_g => setsphericalharmonicbasisgravityconstant, &
                                     sh_clear_basis => clearsphericalharmonicbasis, &
                                     sh_init_basis => initsphericalharmonicbasis, &
                                     sh_default_init_basis => defaultinitsphericalharmonicbasis, &
                                     sh_project_density => project_axisymmetric_density_generic, &
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

    REAL*8, PUBLIC :: GRAVITY_G = GRAVITY_G_DEFAULT
    LOGICAL, PUBLIC :: GRAVITY_G_IS_DEFAULT = .TRUE.
    LOGICAL, PUBLIC :: GRAVITY_FINALIZED = .FALSE.
    INTEGER, PUBLIC :: GRAVITY_NCOMP = 0
    INTEGER, DIMENSION(GRAVITY_MAX_NCOMP), PUBLIC :: GRAVITY_KIND = 0
    REAL*8, DIMENSION(GRAVITY_MAX_PARAMS, GRAVITY_MAX_NCOMP), PUBLIC :: GRAVITY_PARAMS = 0.0D0

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
        GRAVITY_FINALIZED = .TRUE.
    END SUBROUTINE finalizegravity

    SUBROUTINE evaluategravityforces(n, x, y, z, ax, ay, az)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: ax, ay, az
        REAL*8, DIMENSION(n,3) :: force_tmp
        REAL*8, DIMENSION(n) :: ax_c, ay_c, az_c
        INTEGER :: i

        ax = 0.0D0
        ay = 0.0D0
        az = 0.0D0

        IF (.NOT. GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: evaluategravityforces: call finalizegravity first"
            RETURN
        END IF

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
                IF (.NOT. BASIS_EXPANSION_INITIALIZED) THEN
                    IF (.NOT. BASIS_GRID_SET) CALL sh_default_init_basis()
                    CALL sh_project_density(GRAVITY_PARAMS(1:3, i), exponentialoblatehalo_density)
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
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: phi
        REAL*8, DIMENSION(n) :: phi_c
        INTEGER :: i

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
            CASE (GRAVITY_KIND_EXPONENTIALOBLATEHALO)
                IF (.NOT. BASIS_EXPANSION_INITIALIZED) THEN
                    IF (.NOT. BASIS_GRID_SET) CALL sh_default_init_basis()
                    CALL sh_project_density(GRAVITY_PARAMS(1:3, i), exponentialoblatehalo_density)
                    CALL sh_compute_phi_tables()
                END IF
                CALL sh_eval_potential(n, x, y, z, phi_c)
                phi = phi + phi_c
            END SELECT
        END DO
    END SUBROUTINE evaluategravitypotential

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

END MODULE gravity
