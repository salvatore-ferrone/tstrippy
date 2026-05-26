MODULE hostcluster
    ! Contract (table-first, extensible for laws):
    ! 1) configure_hostcluster_kinematics(...) is required
    ! 2) configure_hostcluster_structure(model_name, params) sets constant baseline params
    ! 3) Per-parameter overrides are optional and replace previous values with a warning:
    !      configure_hostcluster_structure_param_table(param_index, times, values)
    !      configure_hostcluster_structure_param_law(param_index, law_name, law_parameters)
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

        SUBROUTINE force_eval_iface(n,x,y,z,force)
            INTEGER, INTENT(IN) :: n
            REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
            REAL*8, INTENT(OUT), DIMENSION(n,3) :: force            
        END SUBROUTINE force_eval_iface

        SUBROUTINE potential_eval_iface(n, x, y, z, phi)
            INTEGER, INTENT(IN) :: n
            REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
            REAL*8, INTENT(OUT), DIMENSION(n) :: phi
        END SUBROUTINE potential_eval_iface


    END INTERFACE

    ! make a derived type to handle the host parameters
    TYPE, PRIVATE :: structural_parameter_t
        CHARACTER(LEN=32) :: name 
         ! 0: constant, 1: table, 2: law
        INTEGER :: evolution_type = 0 ! default at constant 
        REAL*8 :: initial_value ! the value extracted when `configure_hostcluster_structure` is called    
        ! table data
        REAL*8, ALLOCATABLE :: timestamps(:)
        INTEGER :: closest_timestamp = 0 
        REAL*8, ALLOCATABLE :: values(:)
        ! for a law
        ! PROCEDURE(parameter_law_iface), POINTER, NOPASS :: law_eval => NULL() ! default to null
        contains 
            PROCEDURE :: value_at
    END TYPE structural_parameter_t


    ! set the procedure for setting the force and potential evaluator
    PROCEDURE(force_eval_iface), pointer, private :: model_force => NULL()
    PROCEDURE(potential_eval_iface), pointer, private :: model_potential => NULL()

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
    PUBLIC :: configure_hostcluster_structure
    PUBLIC :: finalize_hostcluster
    PUBLIC :: update_hostcluster_state
    PUBLIC :: eval_force
    PUBLIC :: get_kinematics


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
    SUBROUTINE configure_hostcluster_structure(model_name, params, nparams)
        CHARACTER(LEN=*), INTENT(IN) :: model_name
        INTEGER, INTENT(IN) :: nparams
        REAL*8, INTENT(IN), DIMENSION(nparams) :: params

        IF (.NOT. HOST_REGISTERED) THEN
            PRINT*, "WARNING: hostcluster not registered. Call add_hostcluster first"
            RETURN
        END IF

        IF (nparams < 1) THEN
            PRINT*, "WARNING: configure_hostcluster_structure requires nparams >= 1"
            RETURN
        END IF

        IF (HOST_MODEL_SET) THEN
            PRINT*, "WARNING: hostcluster model updated; clearing previous parameter overrides"
        END IF


        HOST_BACKEND_NAME = TRIM(model_name)
        host_params_constant = params
        host_params_current = params

        SELECT CASE (TRIM(MODEL_NAME))
        CASE ("plummer")
            model_force => plummer_force
            model_potential => plummer_potential
        CASE DEFAULT
            PRINT*, "ERROR: unknown hostcluster model: ", TRIM(model_name)
            NULLIFY(model_force)
            NULLIFY(model_potential)
        END SELECT

        HOST_MODEL_SET = .TRUE.
        HOST_NPARAMS = nparams
        HOST_FINALIZED = .FALSE.
    END SUBROUTINE configure_hostcluster_structure

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

    SUBROUTINE update_hostcluster_state(t)
        REAL*8, INTENT(IN) :: t

        IF (.NOT. HOST_FINALIZED) RETURN
        IF (.NOT. ALLOCATED(host_times)) RETURN

        CALL query_current_kinematics(t)

    END SUBROUTINE update_hostcluster_state

    SUBROUTINE eval_force(nparticles, x, y, z, ax, ay, az)
        INTEGER, INTENT(IN) :: nparticles
        REAL*8, INTENT(IN), DIMENSION(nparticles) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(nparticles) :: ax, ay, az
        REAL*8, DIMENSION(nparticles,3) :: force_temp

        REAL*8, DIMENSION(nparticles) :: dx,dy,dz

        dx = x - host_x_current
        dy = y - host_y_current
        dz = z - host_z_current

        call model_force(nparticles,dx,dy,dz,force_temp)
        ax = force_temp(:,1)
        ay = force_temp(:,2)
        az = force_temp(:,3)

    END SUBROUTINE eval_force


    !! TO INTERFACE WITH SIMULATOR
    ! in hostcluster.f90 (inside CONTAINS)
    SUBROUTINE get_kinematics(ntimes, t, x, y, z, vx, vy, vz, ok)
        INTEGER, INTENT(IN) :: ntimes
        REAL*8, INTENT(OUT), DIMENSION(ntimes) :: t, x, y, z, vx, vy, vz
        LOGICAL, INTENT(OUT) :: ok

        ok = .FALSE.

        IF (.NOT. HOST_REGISTERED) THEN
            PRINT*, "WARNING: get_kinematics: host is not registered"
            RETURN
        END IF

        IF (.NOT. ALLOCATED(host_times)) THEN
            PRINT*, "WARNING: get_kinematics: host kinematics are not set"
            RETURN
        END IF

        IF (SIZE(host_times) /= ntimes) THEN
            PRINT*, "WARNING: get_kinematics: ntimes mismatch"
            RETURN
        END IF

        t  = host_times
        x  = host_x
        y  = host_y
        z  = host_z
        vx = host_vx
        vy = host_vy
        vz = host_vz
        ok = .TRUE.
    END SUBROUTINE GET_KINEMATICS    

    !!! MODELS 

    ! ANALYTICAL MODELS
    SUBROUTINE plummer_force(n, x, y, z, force)
        
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n,3) :: force
        REAL*8, DIMENSION(n) :: r, amod
        REAL*8 :: m, b

        m = host_params_current(1)
        b = host_params_current(2)
        r = SQRT(x*x + y*y + z*z)
        amod = -G_hostcluster*m / (r*r + b*b)**1.5

        force(:,1) = amod*x
        force(:,2) = amod*y
        force(:,3) = amod*z
    END SUBROUTINE plummer_force

    SUBROUTINE plummer_potential(n, x, y, z, phi)
        
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: phi
        REAL*8, DIMENSION(n) :: r
        REAL*8 :: m, b

        m = host_params_current(1)
        b = host_params_current(2)
        r = SQRT(x*x + y*y + z*z)
        phi = -G_hostcluster*m / SQRT(r*r + b*b)
    END SUBROUTINE plummer_potential

END MODULE hostcluster
