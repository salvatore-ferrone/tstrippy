MODULE gravity
    USE sphericalharmonicsbfe, ONLY: BASIS_GRID_SET, BASIS_EXPANSION_INITIALIZED, &
                                     setsphericalharmonicbasisgravityconstant, &
                                     clearsphericalharmonicbasis, &
                                     initsphericalharmonicbasis, &
                                     defaultinitsphericalharmonicbasis, &
                                     initsphericalharmoniccomponentphi, &
                                     storesphericalharmoniccomponentphi, &
                                     loadsphericalharmoniccomponentphi, &
                                     project_axisym_density_generic, &
                                     compute_phi_tables_from_rho, &
                                     sphericalharmonicbasisforce, &
                                     sphericalharmonicbasispotential
    USE besselbfe,  bessel_clear                        => clear, &
                    backend_bessel_initialize           => initialize, &
                    bessel_force                        => force, &
                    bessel_potential                    => potential, &
                    bessel_set_gravitational_constant   => set_gravitational_constant, &
                    bessel_allocate_component_tables    => allocate_component_tables, &
                    backend_bessel_set_component_scales => set_component_scales, &
                    bessel_project_density              => project_density, &
                    bessel_default_initialize           => default_initialize, &
                    bessel_load_component               => load_component
    
    IMPLICIT NONE
    
    REAL*8, PARAMETER, PUBLIC :: GRAVITY_G_DEFAULT = 4.30091727D-6
    INTEGER, PARAMETER, PUBLIC :: GRAVITY_MAX_NCOMP = 16
    INTEGER, PARAMETER, PUBLIC :: GRAVITY_MAX_PARAMS = 16

    INTEGER, PARAMETER, PRIVATE :: BACKEND_ANALYTIC = 1
    INTEGER, PARAMETER, PRIVATE :: BACKEND_SH = 2
    INTEGER, PARAMETER, PRIVATE :: BACKEND_BESSEL = 3
    REAL*8, PARAMETER, PRIVATE :: BESSEL_DEFAULT_R_SCALE = 1.0D0
    REAL*8, PARAMETER, PRIVATE :: BESSEL_DEFAULT_Z_SCALE = 1.0D0
    REAL*8, PARAMETER, PRIVATE :: BESSEL_DEFAULT_SCALE_DIVISOR = 3.0D0

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
    CHARACTER(LEN=32), DIMENSION(GRAVITY_MAX_NCOMP), PUBLIC :: COMPONENT_MODEL_NAMES = ''
    LOGICAL, DIMENSION(GRAVITY_MAX_NCOMP), PRIVATE :: BESSEL_SCALE_OVERRIDE_SET = .FALSE.
    REAL*8, DIMENSION(GRAVITY_MAX_NCOMP), PRIVATE :: BESSEL_SCALE_OVERRIDE_R = 0.0D0
    REAL*8, DIMENSION(GRAVITY_MAX_NCOMP), PRIVATE :: BESSEL_SCALE_OVERRIDE_Z = 0.0D0

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
    PRIVATE :: count_sh_components, sh_slot_for_component, count_bessel_components, bessel_slot_for_component

CONTAINS

    SUBROUTINE ensure_component_handlers_initialized()


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
        CALL register_handler_bessel("exponentialdisk", 3, exponentialdisk_density)

        COMPONENT_HANDLERS_INITIALIZED = .TRUE.
    END SUBROUTINE ensure_component_handlers_initialized

    !! SUBROUINES FOR ORGANIZING THE FORCE POINTERS 
    SUBROUTINE register_handler_analytic(model_name, nparams, force_proc, potential_proc)
        
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
        
        INTEGER, INTENT(IN) :: i_comp
        INTEGER :: i_handler

        component_is_sh = .FALSE.
        IF (i_comp < 1 .OR. i_comp > GRAVITY_NCOMP) RETURN
        i_handler = COMPONENT_HANDLER_SLOT(i_comp)
        IF (i_handler < 1 .OR. i_handler > N_COMPONENT_HANDLERS) RETURN

        component_is_sh = (COMPONENT_HANDLERS(i_handler)%backend == BACKEND_SH)
    END FUNCTION component_is_sh

    LOGICAL FUNCTION component_is_bessel(i_comp)
        
        INTEGER, INTENT(IN) :: i_comp
        INTEGER :: i_handler

        component_is_bessel = .FALSE.
        IF (i_comp < 1 .OR. i_comp > GRAVITY_NCOMP) RETURN
        i_handler = COMPONENT_HANDLER_SLOT(i_comp)
        IF (i_handler < 1 .OR. i_handler > N_COMPONENT_HANDLERS) RETURN

        component_is_bessel = (COMPONENT_HANDLERS(i_handler)%backend == BACKEND_BESSEL)
    END FUNCTION component_is_bessel

    ! MODULE-STATE subroutines
    SUBROUTINE clear()
        
        CALL ensure_component_handlers_initialized()
        GRAVITY_G = GRAVITY_G_DEFAULT
        GRAVITY_G_IS_DEFAULT = .TRUE.
        GRAVITY_FINALIZED = .FALSE.
        GRAVITY_NCOMP = 0
        COMPONENT_HANDLER_SLOT = 0
        COMPONENT_MODEL_NAMES = ''
        GRAVITY_PARAMS = 0.0D0
        BESSEL_SCALE_OVERRIDE_SET = .FALSE.
        BESSEL_SCALE_OVERRIDE_R = 0.0D0
        BESSEL_SCALE_OVERRIDE_Z = 0.0D0
        CALL clearsphericalharmonicbasis()
        CALL bessel_clear()
    END SUBROUTINE clear

    SUBROUTINE set_gravitational_constant(g)
        
        REAL*8, INTENT(IN) :: g
        IF (GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: set_gravitational_constant: cannot change G after finalize"
            RETURN
        END IF
        IF (g <= 0.0D0) THEN
            WRITE(*,'(A)') "WARNING: set_gravitational_constant: G must be positive"
            RETURN
        END IF
        GRAVITY_G = g
        GRAVITY_G_IS_DEFAULT = .FALSE.
        CALL setsphericalharmonicbasisgravityconstant(GRAVITY_G)
        call bessel_set_gravitational_constant(GRAVITY_G)
    END SUBROUTINE set_gravitational_constant

    SUBROUTINE add_component(model_name, params, nparams)
        
        CHARACTER(LEN=*), INTENT(IN) :: model_name
        INTEGER, INTENT(IN) :: nparams
        REAL*8, INTENT(IN), DIMENSION(nparams) :: params
        INTEGER :: i_handler

        IF (GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: add_component: cannot add components after finalize"
            RETURN
        END IF
        IF (GRAVITY_NCOMP >= GRAVITY_MAX_NCOMP) THEN
            WRITE(*,'(A)') "WARNING: add_component: maximum number of components reached"
            RETURN
        END IF

        i_handler = handler_index_from_name(model_name)
        IF (i_handler <= 0) THEN
            WRITE(*,'(A)') "WARNING: add_component: unknown model_name", model_name
            RETURN
        END IF

        IF (nparams /= COMPONENT_HANDLERS(i_handler)%nparams) THEN
            WRITE(*,'(A,I0,A,I0)') "WARNING: add_component: expected", COMPONENT_HANDLERS(i_handler)%nparams, &
                                   " params, got", nparams
            RETURN
        END IF

        GRAVITY_NCOMP = GRAVITY_NCOMP + 1
        COMPONENT_HANDLER_SLOT(GRAVITY_NCOMP) = i_handler
        GRAVITY_PARAMS(1:nparams, GRAVITY_NCOMP) = params(1:nparams)
        COMPONENT_MODEL_NAMES(GRAVITY_NCOMP) = TRIM(model_name)
    END SUBROUTINE add_component

    SUBROUTINE finalize()
        
        INTEGER :: i, n_sh, n_bessel, i_sh, i_bessel, i_sh_slot, i_handler
        REAL*8 :: bessel_r_scale, bessel_z_scale
        IF (GRAVITY_NCOMP < 1) THEN
            WRITE(*,'(A)') "WARNING: finalize: no components registered"
            RETURN
        END IF

        n_sh = 0
        n_bessel = 0
        i_sh = 0
        i_bessel = 0
        DO i = 1, GRAVITY_NCOMP
            i_handler = COMPONENT_HANDLER_SLOT(i)
            IF (i_handler < 1 .OR. i_handler > N_COMPONENT_HANDLERS) THEN
                WRITE(*,'(A)') "WARNING: finalize: component slot is uninitialized"
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
            IF (.NOT. BASIS_GRID_SET) CALL defaultinitsphericalharmonicbasis()
            i_handler = COMPONENT_HANDLER_SLOT(i_sh)
            CALL project_axisym_density_generic(GRAVITY_PARAMS(1:COMPONENT_HANDLERS(i_handler)%nparams, i_sh), &
                                    COMPONENT_HANDLERS(i_handler)%density_proc)
            CALL compute_phi_tables_from_rho()
        ELSE IF (n_sh > 1) THEN
            IF (.NOT. BASIS_GRID_SET) CALL defaultinitsphericalharmonicbasis()
            CALL initsphericalharmoniccomponentphi(n_sh)

            i_sh_slot = 0
            DO i = 1, GRAVITY_NCOMP
                IF (component_is_sh(i)) THEN
                    i_sh_slot = i_sh_slot + 1
                    i_handler = COMPONENT_HANDLER_SLOT(i)
                    CALL project_axisym_density_generic(GRAVITY_PARAMS(1:COMPONENT_HANDLERS(i_handler)%nparams, i), &
                                            COMPONENT_HANDLERS(i_handler)%density_proc)
                    CALL compute_phi_tables_from_rho()
                    CALL storesphericalharmoniccomponentphi(i_sh_slot)
                END IF
            END DO
        END IF

        ! Minimal bessel table build
        IF (n_bessel >= 1) THEN
            IF (.NOT. BESSEL_INITIALIZED) CALL bessel_default_initialize()
            CALL bessel_allocate_component_tables(n_bessel)

            i_bessel = 0
            DO i = 1, GRAVITY_NCOMP
                IF (component_is_bessel(i)) THEN
                    i_bessel = i_bessel + 1
                    i_handler = COMPONENT_HANDLER_SLOT(i)

                    ! Auto-select per-component table scales from model params for
                    ! exponentialdisk-style profiles: params=(Sigma0, hR, hZ).
                    ! ASSUMES THAT R SCALE LENGTH IS THE SECOND PARAMETER 
                    ! ASSUMES THAT Z SCALE LENGTH IS THE THIRD PARAMETER !! 
                    ! RESPECT THE CALLING SEQUENCE !
                    bessel_r_scale = BESSEL_DEFAULT_R_SCALE / BESSEL_DEFAULT_SCALE_DIVISOR
                    bessel_z_scale = BESSEL_DEFAULT_Z_SCALE / BESSEL_DEFAULT_SCALE_DIVISOR
                    IF (COMPONENT_HANDLERS(i_handler)%nparams >= 3) THEN
                        IF (GRAVITY_PARAMS(2, i) > 0.0D0) bessel_r_scale = GRAVITY_PARAMS(2, i) / BESSEL_DEFAULT_SCALE_DIVISOR
                        IF (GRAVITY_PARAMS(3, i) > 0.0D0) bessel_z_scale = GRAVITY_PARAMS(3, i) / BESSEL_DEFAULT_SCALE_DIVISOR
                    END IF

                    IF (BESSEL_SCALE_OVERRIDE_SET(i_bessel)) THEN
                        bessel_r_scale = BESSEL_SCALE_OVERRIDE_R(i_bessel)
                        bessel_z_scale = BESSEL_SCALE_OVERRIDE_Z(i_bessel)
                    END IF
                    CALL backend_bessel_set_component_scales(i_bessel, bessel_r_scale, bessel_z_scale)

                    CALL bessel_project_density(i_bessel, &
                                                COMPONENT_HANDLERS(i_handler)%density_proc, &
                                                GRAVITY_PARAMS(1:COMPONENT_HANDLERS(i_handler)%nparams, i) )
                END IF
            END DO
        END IF

        GRAVITY_FINALIZED = .TRUE.
    END SUBROUTINE finalize

    SUBROUTINE ensure_sh_component_tables_loaded(i_comp, n_sh)
        
        INTEGER, INTENT(IN) :: i_comp, n_sh
        INTEGER :: i_handler

        i_handler = COMPONENT_HANDLER_SLOT(i_comp)
        IF (i_handler <= 0) RETURN

        IF (n_sh > 1) THEN
            CALL loadsphericalharmoniccomponentphi(sh_slot_for_component(i_comp))
            RETURN
        END IF

        IF (.NOT. BASIS_EXPANSION_INITIALIZED) THEN
            IF (.NOT. BASIS_GRID_SET) CALL defaultinitsphericalharmonicbasis()
            CALL project_axisym_density_generic(GRAVITY_PARAMS(1:COMPONENT_HANDLERS(i_handler)%nparams, i_comp), &
                                    COMPONENT_HANDLERS(i_handler)%density_proc)
            CALL compute_phi_tables_from_rho()
        END IF
    END SUBROUTINE ensure_sh_component_tables_loaded

    SUBROUTINE eval_component_force(i_comp, n_sh, n, x, y, z, force_c)
        
        INTEGER, INTENT(IN) :: i_comp, n_sh, n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n,3) :: force_c
        INTEGER :: i_handler, i_bessel_slot

        i_handler = COMPONENT_HANDLER_SLOT(i_comp)
        IF (i_handler <= 0) THEN
            force_c = 0.0D0
            RETURN
        END IF

        IF (COMPONENT_HANDLERS(i_handler)%backend == BACKEND_SH) THEN
            CALL ensure_sh_component_tables_loaded(i_comp, n_sh)
        END IF

        IF (COMPONENT_HANDLERS(i_handler)%backend == BACKEND_BESSEL) THEN
            i_bessel_slot = bessel_slot_for_component(i_comp)
            CALL bessel_load_component(i_bessel_slot)
        END IF

        CALL COMPONENT_HANDLERS(i_handler)%force_proc( &
            GRAVITY_PARAMS(1:COMPONENT_HANDLERS(i_handler)%nparams, i_comp), n, x, y, z, force_c)
    END SUBROUTINE eval_component_force

    SUBROUTINE eval_component_potential(i_comp, n_sh, n, x, y, z, phi_c)
        
        INTEGER, INTENT(IN) :: i_comp, n_sh, n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: phi_c
        INTEGER :: i_handler, i_bessel_slot

        i_handler = COMPONENT_HANDLER_SLOT(i_comp)
        IF (i_handler <= 0) THEN
            phi_c = 0.0D0
            RETURN
        END IF

        IF (COMPONENT_HANDLERS(i_handler)%backend == BACKEND_SH) THEN
            CALL ensure_sh_component_tables_loaded(i_comp, n_sh)
        END IF

        IF (COMPONENT_HANDLERS(i_handler)%backend == BACKEND_BESSEL) THEN
            i_bessel_slot = bessel_slot_for_component(i_comp)
            CALL bessel_load_component(i_bessel_slot)
        END IF

        CALL COMPONENT_HANDLERS(i_handler)%potential_proc( &
            GRAVITY_PARAMS(1:COMPONENT_HANDLERS(i_handler)%nparams, i_comp), n, x, y, z, phi_c)
    END SUBROUTINE eval_component_potential

    SUBROUTINE force_components(n, x, y, z, ax_comp, ay_comp, az_comp)
        
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(16, n) :: ax_comp, ay_comp, az_comp
        REAL*8, DIMENSION(n,3) :: force_tmp
        INTEGER :: i, n_sh

        ax_comp = 0.0D0
        ay_comp = 0.0D0
        az_comp = 0.0D0

        IF (.NOT. GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: evaluategravityforcecomponents: call finalize first"
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
        
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: ax, ay, az
        REAL*8, DIMENSION(n,3) :: force_tmp
        INTEGER :: i, n_sh

        ax = 0.0D0
        ay = 0.0D0
        az = 0.0D0

        IF (.NOT. GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: evaluategravityforces: call finalize first"
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
            WRITE(*,'(A)') "WARNING: evaluategravitypotential: call finalize first"
            RETURN
        END IF

        n_sh = count_sh_components()

        DO i = 1, GRAVITY_NCOMP
            CALL eval_component_potential(i, n_sh, n, x, y, z, phi_c)
            phi = phi + phi_c
        END DO
    END SUBROUTINE potential

    SUBROUTINE potential_components(n, x, y, z, phi_comp)
        
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(16, n) :: phi_comp
        REAL*8, DIMENSION(n) :: phi_c
        INTEGER :: i, n_sh

        phi_comp = 0.0D0

        IF (.NOT. GRAVITY_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: evaluategravitypotentialcomponents: call finalize first"
            RETURN
        END IF

        n_sh = count_sh_components()

        DO i = 1, GRAVITY_NCOMP
            CALL eval_component_potential(i, n_sh, n, x, y, z, phi_c)
            phi_comp(i,:) = phi_c
        END DO
    END SUBROUTINE potential_components

    SUBROUTINE sh_force_from_tables(params, n, x, y, z, force)
        
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n,3) :: force

        CALL sphericalharmonicbasisforce(n, x, y, z, force(:,1), force(:,2), force(:,3))
    END SUBROUTINE sh_force_from_tables

    SUBROUTINE sh_potential_from_tables(params, n, x, y, z, phi)
        
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: phi

        CALL sphericalharmonicbasispotential(n, x, y, z, phi)
    END SUBROUTINE sh_potential_from_tables

    ! Spherical harmonics interfacing subroutines
    INTEGER FUNCTION count_sh_components()
        
        INTEGER :: i

        count_sh_components = 0
        DO i = 1, GRAVITY_NCOMP
            IF (component_is_sh(i)) THEN
                count_sh_components = count_sh_components + 1
            END IF
        END DO
    
    END FUNCTION count_sh_components

    INTEGER FUNCTION sh_slot_for_component(i_comp)
        
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

    SUBROUTINE bessel_initialize(nr, nz, nk_build_in)
        ! Wrapper for setting bessel grid/spectral resolution controls.
        
        INTEGER, INTENT(IN) :: nr, nz, nk_build_in
        CALL backend_bessel_initialize(nr, nz, nk_build_in)
    END SUBROUTINE bessel_initialize

    SUBROUTINE bessel_set_component_scales(component_index, r_scale, z_scale)
        ! Optional user override for per-component bessel domain scales.
        
        INTEGER, INTENT(IN) :: component_index
        REAL*8, INTENT(IN) :: r_scale, z_scale

        IF (component_index < 1 .OR. component_index > GRAVITY_MAX_NCOMP) THEN
            WRITE(*,'(A)') "WARNING: bessel_set_component_scales: invalid component_index"
            RETURN
        END IF
        IF (r_scale <= 0.0D0) THEN
            WRITE(*,'(A)') "WARNING: bessel_set_component_scales: r_scale must be positive"
            RETURN
        END IF
        IF (z_scale <= 0.0D0) THEN
            WRITE(*,'(A)') "WARNING: bessel_set_component_scales: z_scale must be positive"
            RETURN
        END IF

        ! Allow pre-finalize configuration by storing overrides by bessel slot.
        IF (.NOT. GRAVITY_FINALIZED) THEN
            BESSEL_SCALE_OVERRIDE_SET(component_index) = .TRUE.
            BESSEL_SCALE_OVERRIDE_R(component_index) = r_scale
            BESSEL_SCALE_OVERRIDE_Z(component_index) = z_scale
            RETURN
        END IF

        CALL backend_bessel_set_component_scales(component_index, r_scale, z_scale)
    END SUBROUTINE bessel_set_component_scales

    SUBROUTINE bessel_force_wrapper(params, n, x, y, z, force)
        
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n,3) :: force
        REAL*8, DIMENSION(n) :: ax, ay, az

        CALL bessel_force(n, x, y, z, ax, ay, az)
        force(:,1) = ax
        force(:,2) = ay
        force(:,3) = az
    END SUBROUTINE bessel_force_wrapper

    SUBROUTINE bessel_potential_wrapper(params, n, x, y, z, phi)
        
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: phi

        CALL bessel_potential(n, x, y, z, phi)
    END SUBROUTINE bessel_potential_wrapper

    INTEGER FUNCTION count_bessel_components()
        
        INTEGER :: i

        count_bessel_components = 0
        DO i = 1, GRAVITY_NCOMP
            IF (component_is_bessel(i)) THEN
                count_bessel_components = count_bessel_components + 1
            END IF
        END DO
    END FUNCTION count_bessel_components

    INTEGER FUNCTION bessel_slot_for_component(i_comp)
        
        INTEGER, INTENT(IN) :: i_comp
        INTEGER :: i

        bessel_slot_for_component = 0
        IF (i_comp < 1 .OR. i_comp > GRAVITY_NCOMP) RETURN

        DO i = 1, i_comp
            IF (component_is_bessel(i)) THEN
                bessel_slot_for_component = bessel_slot_for_component + 1
            END IF
        END DO
    END FUNCTION bessel_slot_for_component

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! ANALYTICAL POTENTIAL MODELS !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! SPHERES
    SUBROUTINE plummer_force(params, n, x, y, z, force)
        
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

    ! PUBLIC ACCESSORS for Python/IO layer
    SUBROUTINE getcomponentmodelname(i_comp, name)
        INTEGER, INTENT(IN) :: i_comp
        CHARACTER(LEN=32), INTENT(OUT) :: name
        INTEGER :: i_handler
        IF (i_comp < 1 .OR. i_comp > GRAVITY_NCOMP) THEN
            name = ""
            RETURN
        END IF
        i_handler = COMPONENT_HANDLER_SLOT(i_comp)
        name = TRIM(COMPONENT_HANDLERS(i_handler)%model_name)
    END SUBROUTINE getcomponentmodelname

    SUBROUTINE getcomponentnparams(i_comp, nparams_out)
        INTEGER, INTENT(IN) :: i_comp
        INTEGER, INTENT(OUT) :: nparams_out
        IF (i_comp < 1 .OR. i_comp > GRAVITY_NCOMP) THEN
            nparams_out = 0
            RETURN
        END IF
        nparams_out = COMPONENT_HANDLERS(COMPONENT_HANDLER_SLOT(i_comp))%nparams
    END SUBROUTINE getcomponentnparams

END MODULE gravity
