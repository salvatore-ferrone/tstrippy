MODULE hostcluster
    ! Contract (table-first, extensible for laws):
    ! 1) configure_hostcluster_kinematics(...) is required
    ! 2) configure_hostcluster_model(model_name, params) sets constant baseline params
    ! 3) Per-parameter overrides are optional and replace previous values with a warning:
    !      configure_hostcluster_model_param_table(param_index, times, values)
    !      configure_hostcluster_model_param_law(param_index, law_name, law_parameters)
    ! 4) Per-parameter precedence at runtime: table > law(stub) > constant
    ! 5) finalize_hostcluster() validates lifecycle/state only (minimal physics policing)
    
    USE mathutils, only : linear_interp_scalar, is_strictly_monotonic, is_strictly_decreasing
    IMPLICIT NONE

    LOGICAL, PUBLIC :: HOST_REGISTERED = .FALSE.
    LOGICAL, PUBLIC :: HOST_FINALIZED = .FALSE.
    LOGICAL, PUBLIC :: HOST_KINEMATICS_SET = .FALSE.
    LOGICAL, PUBLIC :: HOST_MODEL_SET = .FALSE.
    INTEGER, PUBLIC :: HOST_NPARAMS = 0

    CHARACTER(LEN=64), PUBLIC :: HOST_BACKEND_NAME = ""

    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_times
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_x, host_y, host_z
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_vx, host_vy, host_vz

    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_params_current
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_params_constant

    ABSTRACT INTERFACE 
    ! for the different laws that can be used for the time evolution of the structural parameters
        SUBROUTINE parameter_law_iface(lawparams,time,param)
            REAL*8, INTENT(IN), DIMENSION(:) :: lawparams 
            REAL*8, INTENT(IN) :: time 
            REAL*8, INTENT(OUT) :: param
        END subroutine parameter_law_iface
    END INTERFACE

    ! make a derived type to handle the host parameters
    TYPE, PRIVATE :: structural_parameter_t
        CHARACTER(LEN=32) :: name 
         ! 0: constant, 1: table, 2: law
        INTEGER :: evolution_type = 0 ! default at constant 
        REAL*8 :: initial_value ! the value extracted when `configure_hostcluster_model` is called    
        ! table data
        REAL*8, ALLOCATABLE :: timestamps(:)
        INTEGER :: closest_timestamp = 0 
        REAL*8, ALLOCATABLE :: values(:)
        ! for a law
        ! PROCEDURE(parameter_law_iface), POINTER, NOPASS :: law_eval => NULL() ! default to null
        contains 
            PROCEDURE :: value_at
    END TYPE structural_parameter_t




    LOGICAL, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_param_has_table
    LOGICAL, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_param_has_law
    INTEGER, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_param_table_ntimes
    INTEGER, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_param_law_nparams
    CHARACTER(LEN=64), DIMENSION(:), PUBLIC, ALLOCATABLE :: host_param_law_name

    REAL*8, DIMENSION(:,:), PUBLIC, ALLOCATABLE :: host_param_table_times
    REAL*8, DIMENSION(:,:), PUBLIC, ALLOCATABLE :: host_param_table_values
    REAL*8, DIMENSION(:,:), PUBLIC, ALLOCATABLE :: host_param_law_values

    INTEGER, PRIVATE :: host_param_table_capacity = 0
    INTEGER, PRIVATE :: host_param_law_capacity = 0

    INTEGER, PUBLIC :: host_current_kinematics_time_index = 1 
    REAL*8, PUBLIC :: host_x_current = 0.0D0
    REAL*8, PUBLIC :: host_y_current = 0.0D0
    REAL*8, PUBLIC :: host_z_current = 0.0D0
    REAL*8, PUBLIC :: host_vx_current = 0.0D0
    REAL*8, PUBLIC :: host_vy_current = 0.0D0
    REAL*8, PUBLIC :: host_vz_current = 0.0D0
    LOGICAL, PUBLIC :: kinematics_forward_orbit = .TRUE.

    PUBLIC :: clear
    PUBLIC :: add_hostcluster
    PUBLIC :: configure_hostcluster_kinematics
    PUBLIC :: configure_hostcluster_model
    PUBLIC :: configure_hostcluster_model_param_law
    PUBLIC :: configure_hostcluster_model_param_table
    PUBLIC :: finalize_hostcluster
    PUBLIC :: update_hostcluster_state
    PUBLIC :: force_hostcluster_on_particles

    REAL*8, PARAMETER, PRIVATE :: G_DEFAULT = 4.30091727D-6
    REAL*8, PUBLIC :: G_hostcluster = G_DEFAULT
    LOGICAL, PUBLIC :: G_IS_DEFAULT = .TRUE.




CONTAINS



    !!! GENERAL MODULE ROUTINES 
    SUBROUTINE clear()
        IF (ALLOCATED(host_times)) DEALLOCATE(host_times)
        IF (ALLOCATED(host_x)) DEALLOCATE(host_x)
        IF (ALLOCATED(host_y)) DEALLOCATE(host_y)
        IF (ALLOCATED(host_z)) DEALLOCATE(host_z)
        IF (ALLOCATED(host_vx)) DEALLOCATE(host_vx)
        IF (ALLOCATED(host_vy)) DEALLOCATE(host_vy)
        IF (ALLOCATED(host_vz)) DEALLOCATE(host_vz)

        IF (ALLOCATED(host_params_current)) DEALLOCATE(host_params_current)
        IF (ALLOCATED(host_params_constant)) DEALLOCATE(host_params_constant)
        IF (ALLOCATED(host_param_has_table)) DEALLOCATE(host_param_has_table)
        IF (ALLOCATED(host_param_has_law)) DEALLOCATE(host_param_has_law)
        IF (ALLOCATED(host_param_table_ntimes)) DEALLOCATE(host_param_table_ntimes)
        IF (ALLOCATED(host_param_law_nparams)) DEALLOCATE(host_param_law_nparams)
        IF (ALLOCATED(host_param_law_name)) DEALLOCATE(host_param_law_name)
        IF (ALLOCATED(host_param_table_times)) DEALLOCATE(host_param_table_times)
        IF (ALLOCATED(host_param_table_values)) DEALLOCATE(host_param_table_values)
        IF (ALLOCATED(host_param_law_values)) DEALLOCATE(host_param_law_values)

        HOST_REGISTERED = .FALSE.
        HOST_FINALIZED = .FALSE.
        HOST_KINEMATICS_SET = .FALSE.
        HOST_MODEL_SET = .FALSE.
        G_hostcluster = G_DEFAULT
        G_IS_DEFAULT = .TRUE.
        HOST_NPARAMS = 0
        HOST_BACKEND_NAME = ""
        host_x_current = 0.0D0
        host_y_current = 0.0D0
        host_z_current = 0.0D0
        host_vx_current = 0.0D0
        host_vy_current = 0.0D0
        host_vz_current = 0.0D0
        host_param_table_capacity = 0
        host_param_law_capacity = 0
    END SUBROUTINE clear

    SUBROUTINE set_gravitational_constant(g)

        REAL*8, INTENT(IN) :: g
        
        IF (HOST_FINALIZED) THEN
            WRITE(*,'(A)') "WARNING: set_gravitational_constant: cannot change G after finalize"
            RETURN
        END IF
        IF (g <= 0.0D0) THEN
            WRITE(*,'(A)') "WARNING: set_gravitational_constant: G must be positive"
            RETURN
        END IF
        
        G_hostcluster = g
        G_IS_DEFAULT = .FALSE.
    END SUBROUTINE set_gravitational_constant

    SUBROUTINE add_hostcluster()
        HOST_REGISTERED = .TRUE.
        HOST_FINALIZED = .FALSE.
    END SUBROUTINE add_hostcluster

    SUBROUTINE finalize_hostcluster()
        HOST_FINALIZED = .FALSE.

        IF (.NOT. HOST_REGISTERED) THEN
            PRINT*, "WARNING: finalize_hostcluster called before add_hostcluster"
            RETURN
        END IF

        IF (.NOT. HOST_KINEMATICS_SET) THEN
            PRINT*, "WARNING: finalize_hostcluster requires host kinematics"
            RETURN
        END IF

        IF (LEN_TRIM(HOST_BACKEND_NAME) < 1) THEN
            PRINT*, "WARNING: finalize_hostcluster requires backend selection"
            RETURN
        END IF

        IF (.NOT. HOST_MODEL_SET) THEN
            PRINT*, "WARNING: finalize_hostcluster requires model and constant parameters"
            RETURN
        END IF

        CALL query_current_structural_params(host_times(1))
        HOST_FINALIZED = .TRUE.
    END SUBROUTINE finalize_hostcluster


    !!!! HANDELING THE KINEMATICS
    SUBROUTINE configure_hostcluster_kinematics(ntimes, t, x, y, z, vx, vy, vz)
        INTEGER, INTENT(IN) :: ntimes
        REAL*8, INTENT(IN), DIMENSION(ntimes) :: t, x, y, z, vx, vy, vz

        IF (.NOT. HOST_REGISTERED) THEN
            PRINT*, "WARNING: hostcluster not registered. Call add_hostcluster first"
            RETURN
        END IF

        IF (ntimes < 2) THEN
            PRINT*, "WARNING: configure_hostcluster_kinematics requires ntimes >= 2"
            RETURN
        END IF

        IF (.NOT. is_strictly_monotonic(t)) THEN
            PRINT*, "WARNING: configure_hostcluster_kinematics requires strictly monotonic times"
            RETURN
        END IF

        IF (ALLOCATED(host_times)) DEALLOCATE(host_times)
        IF (ALLOCATED(host_x)) DEALLOCATE(host_x)
        IF (ALLOCATED(host_y)) DEALLOCATE(host_y)
        IF (ALLOCATED(host_z)) DEALLOCATE(host_z)
        IF (ALLOCATED(host_vx)) DEALLOCATE(host_vx)
        IF (ALLOCATED(host_vy)) DEALLOCATE(host_vy)
        IF (ALLOCATED(host_vz)) DEALLOCATE(host_vz)
        
        ALLOCATE(host_times(ntimes), host_x(ntimes), host_y(ntimes), host_z(ntimes))
        ALLOCATE(host_vx(ntimes), host_vy(ntimes), host_vz(ntimes))
        
        host_times = t
        host_x = x
        host_y = y
        host_z = z
        host_vx = vx
        host_vy = vy
        host_vz = vz
        
        host_x_current = x(1)
        host_y_current = y(1)
        host_z_current = z(1)
        host_vx_current = vx(1)
        host_vy_current = vy(1)
        host_vz_current = vz(1)
        
        if (is_strictly_decreasing(host_times)) kinematics_forward_orbit=.FALSE.
        
        HOST_KINEMATICS_SET = .TRUE.
        HOST_FINALIZED = .FALSE.
    END SUBROUTINE configure_hostcluster_kinematics

    SUBROUTINE query_current_kinematics(query_time)
        REAL*8, INTENT(IN) :: query_time
        REAL*8 :: T0, TF
        INTEGER :: n
        REAL*8 :: alpha, dt
        LOGICAL :: bracketted

        n = SIZE(host_times)

        if (kinematics_forward_orbit) then 
            IF (query_time <= host_times(1)) THEN
                host_x_current = host_x(1)
                host_y_current = host_y(1)
                host_z_current = host_z(1)
                host_vx_current = host_vx(1)
                host_vy_current = host_vy(1)
                host_vz_current = host_vz(1)
                RETURN
            END IF
            IF (query_time >= host_times(n)) THEN
                host_x_current = host_x(n)
                host_y_current = host_y(n)
                host_z_current = host_z(n)
                host_vx_current = host_vx(n)
                host_vy_current = host_vy(n)
                host_vz_current = host_vz(n)
                RETURN
            END IF
        ELSE 
            IF (query_time >= host_times(1)) THEN
                host_x_current = host_x(1)
                host_y_current = host_y(1)
                host_z_current = host_z(1)
                host_vx_current = host_vx(1)
                host_vy_current = host_vy(1)
                host_vz_current = host_vz(1)
                RETURN
            END IF
            IF (query_time <= host_times(n)) THEN
                host_x_current = host_x(n)
                host_y_current = host_y(n)
                host_z_current = host_z(n)
                host_vx_current = host_vx(n)
                host_vy_current = host_vy(n)
                host_vz_current = host_vz(n)
                RETURN
            END IF            
        END IF 

        ! we expect that the user will query this function at timestamps that progress monotonically
        ! therefore, we will save the most recent timestamp, and only advance it need be.
        ! so the search time for this algorithm is O(1) to O(ntimestamps). Not bad. 
        ! For the interpolation method, it doesn't matter if dt is positive or negative

        ! check the timestamp time
        T0 = host_times(host_current_kinematics_time_index)
        TF = host_times(host_current_kinematics_time_index + 1)
        dt = TF - T0
        if (dt==0d0) THEN 
            print*, "WARNING in query_current_kinematics. dt=0"
        end if 

        bracketted = (query_time.lt.T0).NEQV.(query_time.lt.TF)
        ! NEQV is the same as XOR and simplifies this expression: 
        ! since we don't know if a<b or b<a
        ! ( (a<x) AND (x<b) ) OR ( (x<b) AND (x<a) )


        do while (.NOT.bracketted)

            if (kinematics_forward_orbit) then 
                if (query_time.gt.TF) host_current_kinematics_time_index = host_current_kinematics_time_index + 1
                if (query_time.lt.T0) host_current_kinematics_time_index = host_current_kinematics_time_index - 1 
            else
                if (query_time.lt.TF) host_current_kinematics_time_index = host_current_kinematics_time_index + 1
                if (query_time.gt.t0) host_current_kinematics_time_index = host_current_kinematics_time_index - 1
            END IF 

            T0 = host_times(host_current_kinematics_time_index)
            TF = host_times(host_current_kinematics_time_index + 1)
            dt = TF - T0
            if (dt==0d0) THEN 
                print*, "WARNING in query_current_kinematics. dt=0"
            end if             
            bracketted = (query_time.lt.T0).NEQV.(query_time.lt.TF)
        END DO 
        alpha = (query_time - T0) / dt
        host_x_current = linear_interp_scalar(host_x(host_current_kinematics_time_index),host_x(host_current_kinematics_time_index+1), alpha )
        host_y_current = linear_interp_scalar(host_y(host_current_kinematics_time_index),host_y(host_current_kinematics_time_index+1), alpha )
        host_z_current = linear_interp_scalar(host_z(host_current_kinematics_time_index),host_z(host_current_kinematics_time_index+1), alpha )
        host_vx_current = linear_interp_scalar(host_vx(host_current_kinematics_time_index),host_vx(host_current_kinematics_time_index+1), alpha )
        host_vy_current = linear_interp_scalar(host_vy(host_current_kinematics_time_index),host_vy(host_current_kinematics_time_index+1), alpha )
        host_vz_current = linear_interp_scalar(host_vz(host_current_kinematics_time_index),host_vz(host_current_kinematics_time_index+1), alpha )
    END SUBROUTINE query_current_kinematics


    !!!! ROUTINES FOR HANDLING CHANGING STRUCTURAL PARAMETERS 
    REAL*8 FUNCTION value_at(self, t)
        CLASS(structural_parameter_t), INTENT(IN) :: self 
        REAL*8, INTENT(IN) :: t 

        SELECT CASE (self%evolution_type)
        CASE (0)
            value_at = self%initial_value
        case(1)
            ! interpolate
        CASE DEFAULT
            value_at = self%initial_value
        END SELECT

    end function value_at   

    SUBROUTINE configure_hostcluster_model(model_name, params, nparams)
        CHARACTER(LEN=*), INTENT(IN) :: model_name
        INTEGER, INTENT(IN) :: nparams
        REAL*8, INTENT(IN), DIMENSION(nparams) :: params

        IF (.NOT. HOST_REGISTERED) THEN
            PRINT*, "WARNING: hostcluster not registered. Call add_hostcluster first"
            RETURN
        END IF

        IF (nparams < 1) THEN
            PRINT*, "WARNING: configure_hostcluster_model requires nparams >= 1"
            RETURN
        END IF

        IF (HOST_MODEL_SET) THEN
            PRINT*, "WARNING: hostcluster model updated; clearing previous parameter overrides"
        END IF

        CALL allocate_param_state(nparams)

        HOST_BACKEND_NAME = TRIM(model_name)
        host_params_constant = params
        host_params_current = params
        host_param_has_table = .FALSE.
        host_param_has_law = .FALSE.
        host_param_table_ntimes = 0
        host_param_law_nparams = 0
        host_param_law_name = ""

        HOST_MODEL_SET = .TRUE.
        HOST_NPARAMS = nparams
        HOST_FINALIZED = .FALSE.
    END SUBROUTINE configure_hostcluster_model

    SUBROUTINE configure_hostcluster_model_param_law(param_index, law_name, law_params, nparams)
        INTEGER, INTENT(IN) :: param_index
        CHARACTER(LEN=*), INTENT(IN) :: law_name
        INTEGER, INTENT(IN) :: nparams
        REAL*8, INTENT(IN), DIMENSION(nparams) :: law_params

        IF (.NOT. HOST_MODEL_SET) THEN
            PRINT*, "WARNING: configure_hostcluster_model must be called before parameter overrides"
            RETURN
        END IF

        IF (param_index < 1 .OR. param_index > HOST_NPARAMS) THEN
            PRINT*, "WARNING: configure_hostcluster_model_param_law invalid param_index"
            RETURN
        END IF

        IF (nparams < 1) THEN
            PRINT*, "WARNING: configure_hostcluster_model_param_law requires nparams >= 1"
            RETURN
        END IF

        IF (host_param_has_law(param_index)) THEN
            PRINT*, "WARNING: parameter in time being updated"
        END IF

        CALL ensure_law_capacity(nparams)
        host_param_law_values(:, param_index) = 0.0D0
        host_param_law_values(1:nparams, param_index) = law_params
        host_param_law_nparams(param_index) = nparams
        host_param_law_name(param_index) = TRIM(law_name)
        host_param_has_law(param_index) = .TRUE.

        ! Law execution is intentionally deferred; current value remains constant unless table override exists.
        HOST_FINALIZED = .FALSE.
    END SUBROUTINE configure_hostcluster_model_param_law

    SUBROUTINE configure_hostcluster_model_param_table(param_index, times, values, ntimes)
        INTEGER, INTENT(IN) :: param_index
        INTEGER, INTENT(IN) :: ntimes
        REAL*8, INTENT(IN), DIMENSION(ntimes) :: times, values

        IF (.NOT. HOST_MODEL_SET) THEN
            PRINT*, "WARNING: configure_hostcluster_model must be called before parameter overrides"
            RETURN
        END IF

        IF (param_index < 1 .OR. param_index > HOST_NPARAMS) THEN
            PRINT*, "WARNING: configure_hostcluster_model_param_table invalid param_index"
            RETURN
        END IF

        IF (ntimes < 2) THEN
            PRINT*, "WARNING: configure_hostcluster_model_param_table requires ntimes >= 2"
            RETURN
        END IF

        IF (.NOT. is_strictly_monotonic(times)) THEN
            PRINT*, "WARNING: configure_hostcluster_model_param_table requires strictly monotonic times"
            RETURN
        END IF

        IF (host_param_has_table(param_index)) THEN
            PRINT*, "WARNING: parameter in time being updated"
        END IF

        CALL ensure_table_capacity(ntimes)
        host_param_table_times(:, param_index) = 0.0D0
        host_param_table_values(:, param_index) = 0.0D0
        host_param_table_times(1:ntimes, param_index) = times
        host_param_table_values(1:ntimes, param_index) = values
        host_param_table_ntimes(param_index) = ntimes
        host_param_has_table(param_index) = .TRUE.

        HOST_FINALIZED = .FALSE.
    END SUBROUTINE configure_hostcluster_model_param_table



    SUBROUTINE update_hostcluster_state(t)
        REAL*8, INTENT(IN) :: t

        IF (.NOT. HOST_FINALIZED) RETURN
        IF (.NOT. ALLOCATED(host_times)) RETURN

        CALL query_current_kinematics(t)
        CALL query_current_structural_params(t)

    END SUBROUTINE update_hostcluster_state

    SUBROUTINE force_hostcluster_on_particles(nparticles, x, y, z, ax, ay, az, phi)
        INTEGER, INTENT(IN) :: nparticles
        REAL*8, INTENT(IN), DIMENSION(nparticles) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(nparticles) :: ax, ay, az, phi

        ! Base contract scaffold: backend physics is added in later phases.
        ax = 0.0D0
        ay = 0.0D0
        az = 0.0D0
        phi = 0.0D0
    END SUBROUTINE force_hostcluster_on_particles


    SUBROUTINE query_current_structural_params(t)
        REAL*8, INTENT(IN) :: t
        INTEGER :: i
        ! this could be optomised with a derived type and pointer registery, but not f

        IF (.NOT. ALLOCATED(host_params_current)) RETURN
        host_params_current = host_params_constant

        DO i = 1, HOST_NPARAMS
            IF (host_param_has_table(i)) THEN
                host_params_current(i) = interp_table_value(i, t)
            ELSE IF (host_param_has_law(i)) THEN
                ! Law support is intentionally deferred; keep baseline for now.
            END IF
        END DO
    END SUBROUTINE query_current_structural_params

    REAL*8 FUNCTION interp_table_value(param_index, t)
        INTEGER, INTENT(IN) :: param_index
        REAL*8, INTENT(IN) :: t
        INTEGER :: ntime, j
        REAL*8 :: t0, t1, alpha

        
        interp_table_value = host_params_constant(param_index)
        ntime = host_param_table_ntimes(param_index)
        IF (ntime < 2) RETURN

        IF (t <= host_param_table_times(1, param_index)) THEN
            interp_table_value = host_param_table_values(1, param_index)
            RETURN
        END IF
        IF (t >= host_param_table_times(ntime, param_index)) THEN
            interp_table_value = host_param_table_values(ntime, param_index)
            RETURN
        END IF

        DO j = 1, ntime - 1
            t0 = host_param_table_times(j, param_index)
            t1 = host_param_table_times(j + 1, param_index)
            IF ((t0 <= t .AND. t <= t1) .OR. (t1 <= t .AND. t <= t0)) THEN
                alpha = (t - t0) / (t1 - t0)
                interp_table_value = (1.0D0 - alpha) * host_param_table_values(j, param_index) + &
                                     alpha * host_param_table_values(j + 1, param_index)
                RETURN
            END IF
        END DO
    END FUNCTION interp_table_value


    SUBROUTINE allocate_param_state(nparams)
        INTEGER, INTENT(IN) :: nparams

        IF (ALLOCATED(host_params_constant)) DEALLOCATE(host_params_constant)
        IF (ALLOCATED(host_params_current)) DEALLOCATE(host_params_current)
        IF (ALLOCATED(host_param_has_table)) DEALLOCATE(host_param_has_table)
        IF (ALLOCATED(host_param_has_law)) DEALLOCATE(host_param_has_law)
        IF (ALLOCATED(host_param_table_ntimes)) DEALLOCATE(host_param_table_ntimes)
        IF (ALLOCATED(host_param_law_nparams)) DEALLOCATE(host_param_law_nparams)
        IF (ALLOCATED(host_param_law_name)) DEALLOCATE(host_param_law_name)

        ALLOCATE(host_params_constant(nparams), host_params_current(nparams))
        ALLOCATE(host_param_has_table(nparams), host_param_has_law(nparams))
        ALLOCATE(host_param_table_ntimes(nparams), host_param_law_nparams(nparams))
        ALLOCATE(host_param_law_name(nparams))

        host_param_has_table = .FALSE.
        host_param_has_law = .FALSE.
        host_param_table_ntimes = 0
        host_param_law_nparams = 0
        host_param_law_name = ""

        IF (ALLOCATED(host_param_table_times)) DEALLOCATE(host_param_table_times)
        IF (ALLOCATED(host_param_table_values)) DEALLOCATE(host_param_table_values)
        IF (ALLOCATED(host_param_law_values)) DEALLOCATE(host_param_law_values)
        host_param_table_capacity = 0
        host_param_law_capacity = 0
    END SUBROUTINE allocate_param_state

    SUBROUTINE ensure_table_capacity(ntimes)
        INTEGER, INTENT(IN) :: ntimes
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: tmp_times, tmp_values

        IF (ntimes <= host_param_table_capacity .AND. ALLOCATED(host_param_table_times)) RETURN

        ALLOCATE(tmp_times(ntimes, HOST_NPARAMS), tmp_values(ntimes, HOST_NPARAMS))
        tmp_times = 0.0D0
        tmp_values = 0.0D0

        IF (ALLOCATED(host_param_table_times)) THEN
            tmp_times(1:host_param_table_capacity, :) = host_param_table_times
            tmp_values(1:host_param_table_capacity, :) = host_param_table_values
            DEALLOCATE(host_param_table_times)
            DEALLOCATE(host_param_table_values)
        END IF

        ALLOCATE(host_param_table_times(ntimes, HOST_NPARAMS), host_param_table_values(ntimes, HOST_NPARAMS))
        host_param_table_times = tmp_times
        host_param_table_values = tmp_values
        DEALLOCATE(tmp_times, tmp_values)
        host_param_table_capacity = ntimes
    END SUBROUTINE ensure_table_capacity

    SUBROUTINE ensure_law_capacity(nparams)
        INTEGER, INTENT(IN) :: nparams
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: tmp_law

        IF (nparams <= host_param_law_capacity .AND. ALLOCATED(host_param_law_values)) RETURN

        ALLOCATE(tmp_law(nparams, HOST_NPARAMS))
        tmp_law = 0.0D0
        IF (ALLOCATED(host_param_law_values)) THEN
            tmp_law(1:host_param_law_capacity, :) = host_param_law_values
            DEALLOCATE(host_param_law_values)
        END IF

        ALLOCATE(host_param_law_values(nparams, HOST_NPARAMS))
        host_param_law_values = tmp_law
        DEALLOCATE(tmp_law)
        host_param_law_capacity = nparams
    END SUBROUTINE ensure_law_capacity

    ! FORCES
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
        amod = -G_hostcluster*m / (r*r + b*b)**1.5

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
        phi = -G_hostcluster*m / SQRT(r*r + b*b)
    END SUBROUTINE plummer_potential

END MODULE hostcluster
