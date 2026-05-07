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
    USE besselbfe, ONLY: BESSEL_INITIALIZED, BESSEL_NCOMP, &
                         bessel_set_g => bessel_set_gravity_constant, &
                         bessel_default_init => bessel_default_init, &
                         bessel_init_comp => bessel_init_component_tables, &
                         bessel_project_density => bessel_project_axisym_density_generic, &
                         bessel_load_comp => bessel_load_component, &
                         bessel_eval_force => bessel_eval_force, &
                         bessel_eval_potential => bessel_eval_potential
    IMPLICIT NONE

    REAL*8, PARAMETER, PUBLIC :: GRAVITY_G_DEFAULT = 4.30091727D-6
    INTEGER, PARAMETER, PUBLIC :: GRAVITY_MAX_NCOMP = 16
    INTEGER, PARAMETER, PUBLIC :: GRAVITY_MAX_PARAMS = 16

    INTEGER, PARAMETER, PRIVATE :: BACKEND_ANALYTIC = 1
    INTEGER, PARAMETER, PRIVATE :: BACKEND_SH = 2
    INTEGER, PARAMETER, PRIVATE :: BACKEND_BESSEL = 3

    ABSTRACT INTERFACE
        SUBROUTINE force_eval_iface(params, n, x, y, z, force)
            REAL*8, INTENT(IN), DIMENSION(:) :: params
            INTEGER, INTENT(IN) :: n
            REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
            REAL*8, INTENT(OUT), DIMENSION(n,3) :: force
        END SUBROUTINE force_eval_iface

        SUBROUTINE potential_eval_iface(params, n, x, y, z, phi)
            REAL*8, INTENT(IN), DIMENSION(:) :: params
            INTEGER, INTENT(IN) :: n
            REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
            REAL*8, INTENT(OUT), DIMENSION(n) :: phi
        END SUBROUTINE potential_eval_iface

        SUBROUTINE density_eval_iface(params, n, x, y, z, rho)
            REAL*8, INTENT(IN), DIMENSION(:) :: params
            INTEGER, INTENT(IN) :: n
            REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
            REAL*8, INTENT(OUT), DIMENSION(n) :: rho
        END SUBROUTINE density_eval_iface
    END INTERFACE

    TYPE, PRIVATE :: component_handler_t
        CHARACTER(LEN=32) :: model_name = ""
        INTEGER :: nparams = 0
        INTEGER :: backend = BACKEND_ANALYTIC
        PROCEDURE(force_eval_iface), POINTER, NOPASS :: force_proc => NULL()
        PROCEDURE(potential_eval_iface), POINTER, NOPASS :: potential_proc => NULL()
        PROCEDURE(density_eval_iface), POINTER, NOPASS :: density_proc => NULL()
    END TYPE component_handler_t

    REAL*8, PUBLIC :: GRAVITY_G = GRAVITY_G_DEFAULT
    LOGICAL, PUBLIC :: GRAVITY_G_IS_DEFAULT = .TRUE.
    LOGICAL, PUBLIC :: GRAVITY_FINALIZED = .FALSE.
    INTEGER, PUBLIC :: GRAVITY_NCOMP = 0
    REAL*8, DIMENSION(GRAVITY_MAX_PARAMS, GRAVITY_MAX_NCOMP), PUBLIC :: GRAVITY_PARAMS = 0.0D0
    INTEGER, DIMENSION(GRAVITY_MAX_NCOMP), PRIVATE :: COMPONENT_HANDLER_SLOT = 0

    INTEGER, PARAMETER, PRIVATE :: MAX_COMPONENT_HANDLERS = 16
    TYPE(component_handler_t), DIMENSION(MAX_COMPONENT_HANDLERS), PRIVATE :: COMPONENT_HANDLERS
    INTEGER, PRIVATE :: N_COMPONENT_HANDLERS = 0
    LOGICAL, PRIVATE :: COMPONENT_HANDLERS_INITIALIZED = .FALSE.

    PUBLIC :: ibata2024halo_density
    PRIVATE :: ensure_component_handlers_initialized, register_handler_analytic, register_handler_sh, register_handler_bessel
    PRIVATE :: handler_index_from_name, component_is_sh, component_is_bessel
    PRIVATE :: ensure_sh_component_tables_loaded, eval_component_force, eval_component_potential
    PRIVATE :: sh_force_from_tables, sh_potential_from_tables
    PRIVATE :: bessel_force_wrapper, bessel_potential_wrapper
    PRIVATE :: count_sh_components, sh_slot_for_component

CONTAINS

    SUBROUTINE ensure_component_handlers_initialized()
        IMPLICIT NONE

        IF (COMPONENT_HANDLERS_INITIALIZED) RETURN

        N_COMPONENT_HANDLERS = 0
        CALL register_handler_analytic("plummer", 2, plummer_force, plummer_potential)
        CALL register_handler_analytic("hernquist", 2, hernquist_force, hernquist_potential)
        CALL register_handler_analytic("allensantillianhalo", 4, &
                                       allensantillianhalo_force, allensantillianhalo_potential)
        CALL register_handler_analytic("miyamotonagai", 3, &
                                       miyamotonagai_force, miyamotonagai_potential)
        CALL register_handler_analytic("longmuralibar", 4, &
                                       longmuralibar_force, longmuralibar_potential)
        CALL register_handler_analytic("pouliasis2017pii", 10, &
                                       pouliasis2017pii_force, pouliasis2017pii_potential)
        CALL register_handler_sh("exponentialoblatehalo", 3, &
                                 exponentialoblatehalo_density)
        CALL register_handler_sh("ibata2024halo", 6, ibata2024halo_density)
        CALL register_handler_bessel("exponential_disk_bessel", 2, exponentialdisk_density)

        COMPONENT_HANDLERS_INITIALIZED = .TRUE.
    END SUBROUTINE ensure_component_handlers_initialized

    SUBROUTINE register_handler_analytic(model_name, nparams, force_proc, potential_proc)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: model_name
        INTEGER, INTENT(IN) :: nparams
        PROCEDURE(force_eval_iface) :: force_proc
        PROCEDURE(potential_eval_iface) :: potential_proc

        IF (N_COMPONENT_HANDLERS >= MAX_COMPONENT_HANDLERS) THEN
            WRITE(*,'(A)') "WARNING: register_handler_analytic: exceeded MAX_COMPONENT_HANDLERS"
            RETURN
        END IF

        N_COMPONENT_HANDLERS = N_COMPONENT_HANDLERS + 1
        COMPONENT_HANDLERS(N_COMPONENT_HANDLERS)%model_name = model_name
        COMPONENT_HANDLERS(N_COMPONENT_HANDLERS)%nparams = nparams
        COMPONENT_HANDLERS(N_COMPONENT_HANDLERS)%backend = BACKEND_ANALYTIC
        COMPONENT_HANDLERS(N_COMPONENT_HANDLERS)%force_proc => force_proc
        COMPONENT_HANDLERS(N_COMPONENT_HANDLERS)%potential_proc => potential_proc
    END SUBROUTINE register_handler_analytic

    SUBROUTINE register_handler_sh(model_name, nparams, density_proc)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: model_name
        INTEGER, INTENT(IN) :: nparams
        PROCEDURE(density_eval_iface) :: density_proc

        IF (N_COMPONENT_HANDLERS >= MAX_COMPONENT_HANDLERS) THEN
            WRITE(*,'(A)') "WARNING: register_handler_sh: exceeded MAX_COMPONENT_HANDLERS"
            RETURN
        END IF

        N_COMPONENT_HANDLERS = N_COMPONENT_HANDLERS + 1
        COMPONENT_HANDLERS(N_COMPONENT_HANDLERS)%model_name = model_name
        COMPONENT_HANDLERS(N_COMPONENT_HANDLERS)%nparams = nparams
        COMPONENT_HANDLERS(N_COMPONENT_HANDLERS)%backend = BACKEND_SH
        COMPONENT_HANDLERS(N_COMPONENT_HANDLERS)%force_proc => sh_force_from_tables
        COMPONENT_HANDLERS(N_COMPONENT_HANDLERS)%potential_proc => sh_potential_from_tables
        COMPONENT_HANDLERS(N_COMPONENT_HANDLERS)%density_proc => density_proc
    END SUBROUTINE register_handler_sh

    SUBROUTINE register_handler_bessel(model_name, nparams, density_proc)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: model_name
        INTEGER, INTENT(IN) :: nparams
        PROCEDURE(density_eval_iface) :: density_proc

        IF (N_COMPONENT_HANDLERS >= MAX_COMPONENT_HANDLERS) THEN
            WRITE(*,'(A)') "WARNING: register_handler_bessel: exceeded MAX_COMPONENT_HANDLERS"
            RETURN
        END IF

        N_COMPONENT_HANDLERS = N_COMPONENT_HANDLERS + 1
        COMPONENT_HANDLERS(N_COMPONENT_HANDLERS)%model_name = model_name
        COMPONENT_HANDLERS(N_COMPONENT_HANDLERS)%nparams = nparams
        COMPONENT_HANDLERS(N_COMPONENT_HANDLERS)%backend = BACKEND_BESSEL
        COMPONENT_HANDLERS(N_COMPONENT_HANDLERS)%force_proc => bessel_force_wrapper
        COMPONENT_HANDLERS(N_COMPONENT_HANDLERS)%potential_proc => bessel_potential_wrapper
        COMPONENT_HANDLERS(N_COMPONENT_HANDLERS)%density_proc => density_proc
    END SUBROUTINE register_handler_bessel

    INTEGER FUNCTION handler_index_from_name(model_name)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: model_name
        INTEGER :: i

        CALL ensure_component_handlers_initialized()
        handler_index_from_name = 0
        DO i = 1, N_COMPONENT_HANDLERS
            IF (TRIM(model_name) == TRIM(COMPONENT_HANDLERS(i)%model_name)) THEN
                handler_index_from_name = i
                RETURN
            END IF
        END DO
    END FUNCTION handler_index_from_name

    LOGICAL FUNCTION component_is_sh(i_comp)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: i_comp
        INTEGER :: i_handler

        component_is_sh = .FALSE.
        IF (i_comp < 1 .OR. i_comp > GRAVITY_NCOMP) RETURN
        i_handler = COMPONENT_HANDLER_SLOT(i_comp)
        IF (i_handler < 1 .OR. i_handler > N_COMPONENT_HANDLERS) RETURN

        component_is_sh = (COMPONENT_HANDLERS(i_handler)%backend == BACKEND_SH)
    END FUNCTION component_is_sh

    LOGICAL FUNCTION component_is_bessel(i_comp)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: i_comp
        INTEGER :: i_handler

        component_is_bessel = .FALSE.
        IF (i_comp < 1 .OR. i_comp > GRAVITY_NCOMP) RETURN
        i_handler = COMPONENT_HANDLER_SLOT(i_comp)
        IF (i_handler < 1 .OR. i_handler > N_COMPONENT_HANDLERS) RETURN

        component_is_bessel = (COMPONENT_HANDLERS(i_handler)%backend == BACKEND_BESSEL)
    END FUNCTION component_is_bessel

    ! MODULE-STATE subroutines
    SUBROUTINE cleargravity()
        IMPLICIT NONE
        CALL ensure_component_handlers_initialized()
        GRAVITY_G = GRAVITY_G_DEFAULT
        GRAVITY_G_IS_DEFAULT = .TRUE.
        GRAVITY_FINALIZED = .FALSE.
        GRAVITY_NCOMP = 0
        COMPONENT_HANDLER_SLOT = 0
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
        INTEGER :: i_handler

        IF (GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: addgravitycomponent: cannot add components after finalizegravity"
            RETURN
        END IF
        IF (GRAVITY_NCOMP >= GRAVITY_MAX_NCOMP) THEN
            WRITE(*,'(A)') "WARNING: addgravitycomponent: maximum number of components reached"
            RETURN
        END IF

        i_handler = handler_index_from_name(model_name)
        IF (i_handler <= 0) THEN
            WRITE(*,'(A)') "WARNING: addgravitycomponent: unknown model_name", model_name
            RETURN
        END IF

        IF (nparams /= COMPONENT_HANDLERS(i_handler)%nparams) THEN
            WRITE(*,'(A,I0,A,I0)') "WARNING: addgravitycomponent: expected", COMPONENT_HANDLERS(i_handler)%nparams, &
                                   " params, got", nparams
            RETURN
        END IF

        GRAVITY_NCOMP = GRAVITY_NCOMP + 1
        COMPONENT_HANDLER_SLOT(GRAVITY_NCOMP) = i_handler
        GRAVITY_PARAMS(1:nparams, GRAVITY_NCOMP) = params(1:nparams)
    END SUBROUTINE addgravitycomponent

    SUBROUTINE finalizegravity()
        IMPLICIT NONE
        INTEGER :: i, n_sh, n_bessel, i_sh, i_bessel, i_sh_slot, i_handler
        IF (GRAVITY_NCOMP < 1) THEN
            WRITE(*,'(A)') "WARNING: finalizegravity: no components registered"
            RETURN
        END IF

        n_sh = 0
        n_bessel = 0
        i_sh = 0
        i_bessel = 0
        DO i = 1, GRAVITY_NCOMP
            i_handler = COMPONENT_HANDLER_SLOT(i)
            IF (i_handler < 1 .OR. i_handler > N_COMPONENT_HANDLERS) THEN
                WRITE(*,'(A)') "WARNING: finalizegravity: component slot is uninitialized"
                RETURN
            END IF
            IF (component_is_sh(i)) THEN
                n_sh = n_sh + 1
                i_sh = i
            END IF
            IF (component_is_bessel(i)) THEN
                n_bessel = n_bessel + 1
                i_bessel = i
            END IF
        END DO

        ! Microstep: eager SH table build in finalize for the single-SH-component case.
        IF (n_sh == 1) THEN
            IF (.NOT. BASIS_GRID_SET) CALL sh_default_init_basis()
            i_handler = COMPONENT_HANDLER_SLOT(i_sh)
            CALL sh_project_density(GRAVITY_PARAMS(1:COMPONENT_HANDLERS(i_handler)%nparams, i_sh), &
                                    COMPONENT_HANDLERS(i_handler)%density_proc)
            CALL sh_compute_phi_tables()
        ELSE IF (n_sh > 1) THEN
            IF (.NOT. BASIS_GRID_SET) CALL sh_default_init_basis()
            CALL sh_init_component_phi(n_sh)

            i_sh_slot = 0
            DO i = 1, GRAVITY_NCOMP
                IF (component_is_sh(i)) THEN
                    i_sh_slot = i_sh_slot + 1
                    i_handler = COMPONENT_HANDLER_SLOT(i)
                    CALL sh_project_density(GRAVITY_PARAMS(1:COMPONENT_HANDLERS(i_handler)%nparams, i), &
                                            COMPONENT_HANDLERS(i_handler)%density_proc)
                    CALL sh_compute_phi_tables()
                    CALL sh_store_component_phi(i_sh_slot)
                END IF
            END DO
        END IF

        ! Minimal bessel table build
        IF (n_bessel >= 1) THEN
            IF (.NOT. BESSEL_INITIALIZED) CALL bessel_default_init()
            CALL bessel_init_comp(n_bessel)
            DO i = 1, GRAVITY_NCOMP
                IF (component_is_bessel(i)) THEN
                    i_handler = COMPONENT_HANDLER_SLOT(i)
                    CALL bessel_project_density(i, GRAVITY_PARAMS(1:COMPONENT_HANDLERS(i_handler)%nparams, i), &
                                                COMPONENT_HANDLERS(i_handler)%density_proc)
                END IF
            END DO
        END IF

        GRAVITY_FINALIZED = .TRUE.
    END SUBROUTINE finalizegravity

    SUBROUTINE ensure_sh_component_tables_loaded(i_comp, n_sh)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: i_comp, n_sh
        INTEGER :: i_handler

        i_handler = COMPONENT_HANDLER_SLOT(i_comp)
        IF (i_handler <= 0) RETURN

        IF (n_sh > 1) THEN
            CALL sh_load_component_phi(sh_slot_for_component(i_comp))
            RETURN
        END IF

        IF (.NOT. BASIS_EXPANSION_INITIALIZED) THEN
            IF (.NOT. BASIS_GRID_SET) CALL sh_default_init_basis()
            CALL sh_project_density(GRAVITY_PARAMS(1:COMPONENT_HANDLERS(i_handler)%nparams, i_comp), &
                                    COMPONENT_HANDLERS(i_handler)%density_proc)
            CALL sh_compute_phi_tables()
        END IF
    END SUBROUTINE ensure_sh_component_tables_loaded

    SUBROUTINE eval_component_force(i_comp, n_sh, n, x, y, z, force_c)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: i_comp, n_sh, n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n,3) :: force_c
        INTEGER :: i_handler

        i_handler = COMPONENT_HANDLER_SLOT(i_comp)
        IF (i_handler <= 0) THEN
            force_c = 0.0D0
            RETURN
        END IF

        IF (COMPONENT_HANDLERS(i_handler)%backend == BACKEND_SH) THEN
            CALL ensure_sh_component_tables_loaded(i_comp, n_sh)
        END IF

        IF (COMPONENT_HANDLERS(i_handler)%backend == BACKEND_BESSEL) THEN
            CALL bessel_load_comp(i_comp)
        END IF

        CALL COMPONENT_HANDLERS(i_handler)%force_proc( &
            GRAVITY_PARAMS(1:COMPONENT_HANDLERS(i_handler)%nparams, i_comp), n, x, y, z, force_c)
    END SUBROUTINE eval_component_force

    SUBROUTINE eval_component_potential(i_comp, n_sh, n, x, y, z, phi_c)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: i_comp, n_sh, n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: phi_c
        INTEGER :: i_handler

        i_handler = COMPONENT_HANDLER_SLOT(i_comp)
        IF (i_handler <= 0) THEN
            phi_c = 0.0D0
            RETURN
        END IF

        IF (COMPONENT_HANDLERS(i_handler)%backend == BACKEND_SH) THEN
            CALL ensure_sh_component_tables_loaded(i_comp, n_sh)
        END IF

        IF (COMPONENT_HANDLERS(i_handler)%backend == BACKEND_BESSEL) THEN
            CALL bessel_load_comp(i_comp)
        END IF

        CALL COMPONENT_HANDLERS(i_handler)%potential_proc( &
            GRAVITY_PARAMS(1:COMPONENT_HANDLERS(i_handler)%nparams, i_comp), n, x, y, z, phi_c)
    END SUBROUTINE eval_component_potential

    SUBROUTINE force_components(n, x, y, z, ax_comp, ay_comp, az_comp)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(16, n) :: ax_comp, ay_comp, az_comp
        REAL*8, DIMENSION(n,3) :: force_tmp
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
            CALL eval_component_force(i, n_sh, n, x, y, z, force_tmp)
            ax_comp(i,:) = force_tmp(:,1)
            ay_comp(i,:) = force_tmp(:,2)
            az_comp(i,:) = force_tmp(:,3)
        END DO
    END SUBROUTINE force_components

    SUBROUTINE force(n, x, y, z, ax, ay, az)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: ax, ay, az
        REAL*8, DIMENSION(n,3) :: force_tmp
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
            CALL eval_component_force(i, n_sh, n, x, y, z, force_tmp)
            ax = ax + force_tmp(:,1)
            ay = ay + force_tmp(:,2)
            az = az + force_tmp(:,3)
        END DO
    END SUBROUTINE force

    SUBROUTINE potential(n, x, y, z, phi)
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
            CALL eval_component_potential(i, n_sh, n, x, y, z, phi_c)
            phi = phi + phi_c
        END DO
    END SUBROUTINE potential

    SUBROUTINE potential_components(n, x, y, z, phi_comp)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(16, n) :: phi_comp
        REAL*8, DIMENSION(n) :: phi_c
        INTEGER :: i, n_sh

        phi_comp = 0.0D0

        IF (.NOT. GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: evaluategravitypotentialcomponents: call finalizegravity first"
            RETURN
        END IF

        n_sh = count_sh_components()

        DO i = 1, GRAVITY_NCOMP
            CALL eval_component_potential(i, n_sh, n, x, y, z, phi_c)
            phi_comp(i,:) = phi_c
        END DO
    END SUBROUTINE potential_components

    SUBROUTINE sh_force_from_tables(params, n, x, y, z, force)
        IMPLICIT NONE
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n,3) :: force

        CALL sh_eval_force(n, x, y, z, force(:,1), force(:,2), force(:,3))
    END SUBROUTINE sh_force_from_tables

    SUBROUTINE sh_potential_from_tables(params, n, x, y, z, phi)
        IMPLICIT NONE
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: phi

        CALL sh_eval_potential(n, x, y, z, phi)
    END SUBROUTINE sh_potential_from_tables

    SUBROUTINE bessel_force_wrapper(params, n, x, y, z, force)
        IMPLICIT NONE
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n,3) :: force
        REAL*8, DIMENSION(n) :: ax, ay, az

        CALL bessel_eval_force(n, x, y, z, ax, ay, az)
        force(:,1) = ax
        force(:,2) = ay
        force(:,3) = az
    END SUBROUTINE bessel_force_wrapper

    SUBROUTINE bessel_potential_wrapper(params, n, x, y, z, phi)
        IMPLICIT NONE
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: phi

        CALL bessel_eval_potential(n, x, y, z, phi)
    END SUBROUTINE bessel_potential_wrapper

    ! Spherical harmonics interfacing subroutines
    INTEGER FUNCTION count_sh_components()
        IMPLICIT NONE
        INTEGER :: i

        count_sh_components = 0
        DO i = 1, GRAVITY_NCOMP
            IF (component_is_sh(i)) THEN
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
            IF (component_is_sh(i)) THEN
                sh_slot_for_component = sh_slot_for_component + 1
            END IF
        END DO
    END FUNCTION sh_slot_for_component

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! ANALYTICAL POTENTIAL MODELS !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! SPHERES
    SUBROUTINE plummer_force(params, n, x, y, z, force)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(:) :: params
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
        REAL*8, INTENT(IN), DIMENSION(:) :: params
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
        REAL*8, INTENT(IN), DIMENSION(:) :: params
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
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        REAL*8, INTENT(OUT), DIMENSION(n) :: phi
        REAL*8, DIMENSION(n) :: r
        REAL*8 :: m, a

        m = params(1)
        a = params(2)
        r = SQRT(x*x + y*y + z*z)
        phi = -GRAVITY_G*m / (r + a)
    END SUBROUTINE hernquist_potential

    SUBROUTINE allensantillianhalo_force(params, N, x, y, z, force)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(:) :: params
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
    END SUBROUTINE allensantillianhalo_force

    SUBROUTINE allensantillianhalo_potential(params, N, x, y, z, phi)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(:) :: params
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
    END SUBROUTINE allensantillianhalo_potential

    ! DISKS
    SUBROUTINE miyamotonagai_force(params, N, x, y, z, force)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(:) :: params
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
    END SUBROUTINE miyamotonagai_force

    SUBROUTINE miyamotonagai_potential(params, N, x, y, z, phi)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        REAL*8, INTENT(OUT), DIMENSION(N) :: phi
        REAL*8, DIMENSION(N) :: R, zmod
        REAL*8 :: M, a, b

        M = params(1)
        a = params(2)
        b = params(3)
        R = SQRT(x*x + y*y)
        zmod = a + SQRT(z*z + b*b)
        phi = -GRAVITY_G*M / SQRT(R*R + zmod*zmod)
    END SUBROUTINE miyamotonagai_potential

    ! TRIAXIAL 
    SUBROUTINE longmuralibar_force(params, N, x, y, z, force)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(:) :: params
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
    END SUBROUTINE longmuralibar_force

    SUBROUTINE longmuralibar_potential(params, N, x, y, z, phi)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(:) :: params
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
    END SUBROUTINE longmuralibar_potential
    
    ! COMPOSITE
    SUBROUTINE pouliasis2017pii_force(params, N, x, y, z, force)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        REAL*8, INTENT(OUT), DIMENSION(N,3) :: force
        REAL*8, DIMENSION(3) :: thindisk, thickdisk
        REAL*8, DIMENSION(4) :: halo
        REAL*8, DIMENSION(N,3) :: force_h, force_d1, force_d2

        halo = (/params(1), params(2), params(3), params(4)/)
        thindisk = (/params(5), params(6), params(7)/)
        thickdisk = (/params(8), params(9), params(10)/)

        CALL allensantillianhalo_force(halo, N, x, y, z, force_h)
        CALL miyamotonagai_force(thindisk, N, x, y, z, force_d1)
        CALL miyamotonagai_force(thickdisk, N, x, y, z, force_d2)
        force = force_h + force_d1 + force_d2
    END SUBROUTINE pouliasis2017pii_force

    SUBROUTINE pouliasis2017pii_potential(params, N, x, y, z, phi)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        REAL*8, INTENT(OUT), DIMENSION(N) :: phi
        REAL*8, DIMENSION(3) :: thindisk, thickdisk
        REAL*8, DIMENSION(4) :: halo
        REAL*8, DIMENSION(N) :: phi_h, phi_d1, phi_d2

        halo = (/params(1), params(2), params(3), params(4)/)
        thindisk = (/params(5), params(6), params(7)/)
        thickdisk = (/params(8), params(9), params(10)/)

        CALL allensantillianhalo_potential(halo, N, x, y, z, phi_h)
        CALL miyamotonagai_potential(thindisk, N, x, y, z, phi_d1)
        CALL miyamotonagai_potential(thickdisk, N, x, y, z, phi_d2)
        phi = phi_h + phi_d1 + phi_d2
    END SUBROUTINE pouliasis2017pii_potential


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! DENSITY ONLY PROFILES !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!! AXIS SYMMETRIC !!!
    ! HALOS (spherical harmonics)
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

    ! DISKS (bessel functions)
    SUBROUTINE exponentialdisk_density(params, n, x, y, z, rho)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: rho
        REAL*8 :: sigma0, hR, hZ, rho0
        REAL*8, DIMENSION(n) :: R 

        sigma0  = params(1)
        hR      = params(2)
        hZ      = params(3)
        rho0    = sigma0 / (2 * hZ)
        R = sqrt(x**2 + y**2)

        rho = rho0 * exp( -(R/hR) - abs(z)/hZ)

    end subroutine exponentialdisk_density     

END MODULE gravity
