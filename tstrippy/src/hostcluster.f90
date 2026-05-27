MODULE hostcluster
    ! Contract (table-first, extensible for laws):
    ! 1) configure_hostcluster_kinematics(...) is required
    ! 2) configure_hostcluster_structure(model_name, params) sets constant baseline params
    ! 3) Per-parameter overrides are optional and replace previous values with a warning:
    !      configure_hostcluster_structure_param_table(param_index, times, values)
    !      configure_hostcluster_structure_param_law(param_index, law_name, law_parameters)
    ! 4) Per-parameter precedence at runtime: table > law(stub) > constant
    ! 5) finalize_hostcluster() validates lifecycle/state only (minimal physics policing)
    
    USE mathutils, only : linear_interp_scalar, is_strictly_monotonic, is_strictly_decreasing, bracketed_index_search
    IMPLICIT NONE

    ! MODULE STATE VARIABLES 
    LOGICAL, PUBLIC :: HOST_REGISTERED = .FALSE.
    LOGICAL, PUBLIC :: HOST_FINALIZED = .FALSE.
    LOGICAL, PUBLIC :: HOST_KINEMATICS_SET = .FALSE.
    LOGICAL, PUBLIC :: HOST_STRUCTURE_SET = .FALSE.
    LOGICAL, PUBLIC :: KINEMATICS_FORWARD_ORBIT = .TRUE.

    REAL*8, PARAMETER, PRIVATE :: G_DEFAULT = 4.30091727D-6
    REAL*8, PUBLIC  :: G_HOSTCLUSTER = G_DEFAULT
    LOGICAL, PUBLIC :: G_IS_DEFAULT = .TRUE.

    ! for the different laws that can be used for the time evolution of the structural parameters
    ABSTRACT INTERFACE 
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
        CHARACTER(LEN=64) :: name 
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

    ! FOR THE STRUCTURAL PARAMETER DEFAULTS
    INTEGER, PUBLIC :: HOST_NPARAMS = 0
    CHARACTER(LEN=64), PUBLIC :: HOST_MODEL_NAME = ""
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: HOST_PARAMS_CURRENT
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: HOST_PARAMS_CONSTANT        


    ! THE VARIABLES FOR THE KINEMATICS
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: HOST_TIMES
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: HOST_X, HOST_Y, HOST_Z
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: HOST_VX, HOST_VY, HOST_VZ

    INTEGER, PUBLIC :: HOST_CURRENT_KINEMATICS_TIME_INDEX = 1 
    REAL*8, PUBLIC  :: HOST_X_CURRENT = 0.0D0
    REAL*8, PUBLIC  :: HOST_Y_CURRENT = 0.0D0
    REAL*8, PUBLIC  :: HOST_Z_CURRENT = 0.0D0
    REAL*8, PUBLIC  :: HOST_VX_CURRENT = 0.0D0
    REAL*8, PUBLIC  :: HOST_VY_CURRENT = 0.0D0
    REAL*8, PUBLIC  :: HOST_VZ_CURRENT = 0.0D0

    PUBLIC :: clear
    PUBLIC :: add_hostcluster
    PUBLIC :: configure_hostcluster_kinematics
    PUBLIC :: configure_hostcluster_structure
    PUBLIC :: finalize_hostcluster
    PUBLIC :: update_hostcluster_state
    PUBLIC :: eval_force
    PUBLIC :: get_kinematics


CONTAINS

    !!! GENERAL MODULE ROUTINES 
    SUBROUTINE clear()
        IF (ALLOCATED(HOST_TIMES)) DEALLOCATE(HOST_TIMES)
        IF (ALLOCATED(HOST_X)) DEALLOCATE(HOST_X)
        IF (ALLOCATED(HOST_Y)) DEALLOCATE(HOST_Y)
        IF (ALLOCATED(HOST_Z)) DEALLOCATE(HOST_Z)
        IF (ALLOCATED(HOST_VX)) DEALLOCATE(HOST_VX)
        IF (ALLOCATED(HOST_VY)) DEALLOCATE(HOST_VY)
        IF (ALLOCATED(HOST_VZ)) DEALLOCATE(HOST_VZ)

        IF (ALLOCATED(HOST_PARAMS_CURRENT)) DEALLOCATE(HOST_PARAMS_CURRENT)
        IF (ALLOCATED(HOST_PARAMS_CONSTANT)) DEALLOCATE(HOST_PARAMS_CONSTANT)


        HOST_REGISTERED = .FALSE.
        HOST_FINALIZED = .FALSE.
        HOST_KINEMATICS_SET = .FALSE.
        HOST_STRUCTURE_SET = .FALSE.
        G_HOSTCLUSTER = G_DEFAULT
        G_IS_DEFAULT = .TRUE.
        HOST_NPARAMS = 0
        HOST_MODEL_NAME = ""
        HOST_X_CURRENT = 0.0D0
        HOST_Y_CURRENT = 0.0D0
        HOST_Z_CURRENT = 0.0D0
        HOST_VX_CURRENT = 0.0D0
        HOST_VY_CURRENT = 0.0D0
        HOST_VZ_CURRENT = 0.0D0
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
        
        G_HOSTCLUSTER = g
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

        IF (LEN_TRIM(HOST_MODEL_NAME) < 1) THEN
            PRINT*, "WARNING: finalize_hostcluster requires backend selection"
            RETURN
        END IF

        IF (.NOT. HOST_STRUCTURE_SET) THEN
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

        IF (ALLOCATED(HOST_TIMES)) DEALLOCATE(HOST_TIMES)
        IF (ALLOCATED(HOST_X)) DEALLOCATE(HOST_X)
        IF (ALLOCATED(HOST_Y)) DEALLOCATE(HOST_Y)
        IF (ALLOCATED(HOST_Z)) DEALLOCATE(HOST_Z)
        IF (ALLOCATED(HOST_VX)) DEALLOCATE(HOST_VX)
        IF (ALLOCATED(HOST_VY)) DEALLOCATE(HOST_VY)
        IF (ALLOCATED(HOST_VZ)) DEALLOCATE(HOST_VZ)
        
        ALLOCATE(HOST_TIMES(ntimes), HOST_X(ntimes), HOST_Y(ntimes), HOST_Z(ntimes))
        ALLOCATE(HOST_VX(ntimes), HOST_VY(ntimes), HOST_VZ(ntimes))
        
        HOST_TIMES = t
        HOST_X = x
        HOST_Y = y
        HOST_Z = z
        HOST_VX = vx
        HOST_VY = vy
        HOST_VZ = vz
        
        HOST_X_CURRENT = x(1)
        HOST_Y_CURRENT = y(1)
        HOST_Z_CURRENT = z(1)
        HOST_VX_CURRENT = vx(1)
        HOST_VY_CURRENT = vy(1)
        HOST_VZ_CURRENT = vz(1)
        
        if (is_strictly_decreasing(HOST_TIMES)) KINEMATICS_FORWARD_ORBIT=.FALSE.
        
        HOST_KINEMATICS_SET = .TRUE.
        HOST_FINALIZED = .FALSE.
    END SUBROUTINE configure_hostcluster_kinematics

    SUBROUTINE query_current_kinematics(query_time)
        REAL*8, INTENT(IN) :: query_time
        REAL*8 :: T0, TF
        INTEGER :: n
        REAL*8 :: alpha, dt
        LOGICAL :: bracketted

        n = SIZE(HOST_TIMES)

        if (KINEMATICS_FORWARD_ORBIT) then 
            IF (query_time <= HOST_TIMES(1)) THEN
                HOST_X_CURRENT = HOST_X(1)
                HOST_Y_CURRENT = HOST_Y(1)
                HOST_Z_CURRENT = HOST_Z(1)
                HOST_VX_CURRENT = HOST_VX(1)
                HOST_VY_CURRENT = HOST_VY(1)
                HOST_VZ_CURRENT = HOST_VZ(1)
                RETURN
            END IF
            IF (query_time >= HOST_TIMES(n)) THEN
                HOST_X_CURRENT = HOST_X(n)
                HOST_Y_CURRENT = HOST_Y(n)
                HOST_Z_CURRENT = HOST_Z(n)
                HOST_VX_CURRENT = HOST_VX(n)
                HOST_VY_CURRENT = HOST_VY(n)
                HOST_VZ_CURRENT = HOST_VZ(n)
                RETURN
            END IF
        ELSE 
            IF (query_time >= HOST_TIMES(1)) THEN
                HOST_X_CURRENT = HOST_X(1)
                HOST_Y_CURRENT = HOST_Y(1)
                HOST_Z_CURRENT = HOST_Z(1)
                HOST_VX_CURRENT = HOST_VX(1)
                HOST_VY_CURRENT = HOST_VY(1)
                HOST_VZ_CURRENT = HOST_VZ(1)
                RETURN
            END IF
            IF (query_time <= HOST_TIMES(n)) THEN
                HOST_X_CURRENT = HOST_X(n)
                HOST_Y_CURRENT = HOST_Y(n)
                HOST_Z_CURRENT = HOST_Z(n)
                HOST_VX_CURRENT = HOST_VX(n)
                HOST_VY_CURRENT = HOST_VY(n)
                HOST_VZ_CURRENT = HOST_VZ(n)
                RETURN
            END IF            
        END IF 

        ! do quick search, which which will be between O(1) to O(N_TIME_STAMPS), works if forward or backward
        HOST_CURRENT_KINEMATICS_TIME_INDEX = bracketed_index_search(query_time, HOST_CURRENT_KINEMATICS_TIME_INDEX, HOST_TIMES)

        alpha = (query_time - T0) / dt
        HOST_X_CURRENT = linear_interp_scalar(HOST_X(HOST_CURRENT_KINEMATICS_TIME_INDEX),HOST_X(HOST_CURRENT_KINEMATICS_TIME_INDEX+1), alpha )
        HOST_Y_CURRENT = linear_interp_scalar(HOST_Y(HOST_CURRENT_KINEMATICS_TIME_INDEX),HOST_Y(HOST_CURRENT_KINEMATICS_TIME_INDEX+1), alpha )
        HOST_Z_CURRENT = linear_interp_scalar(HOST_Z(HOST_CURRENT_KINEMATICS_TIME_INDEX),HOST_Z(HOST_CURRENT_KINEMATICS_TIME_INDEX+1), alpha )
        HOST_VX_CURRENT = linear_interp_scalar(HOST_VX(HOST_CURRENT_KINEMATICS_TIME_INDEX),HOST_VX(HOST_CURRENT_KINEMATICS_TIME_INDEX+1), alpha )
        HOST_VY_CURRENT = linear_interp_scalar(HOST_VY(HOST_CURRENT_KINEMATICS_TIME_INDEX),HOST_VY(HOST_CURRENT_KINEMATICS_TIME_INDEX+1), alpha )
        HOST_VZ_CURRENT = linear_interp_scalar(HOST_VZ(HOST_CURRENT_KINEMATICS_TIME_INDEX),HOST_VZ(HOST_CURRENT_KINEMATICS_TIME_INDEX+1), alpha )
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

        IF (HOST_STRUCTURE_SET) THEN
            PRINT*, "WARNING: hostcluster model updated; clearing previous parameter overrides"
        END IF

        SELECT CASE (TRIM(MODEL_NAME))
        CASE ("plummer")
            model_force => plummer_force
            model_potential => plummer_potential
        CASE DEFAULT
            PRINT*, "ERROR: unknown hostcluster model: ", TRIM(model_name)
            NULLIFY(model_force)
            NULLIFY(model_potential)
        END SELECT

        if (ALLOCATED(HOST_PARAMS_CONSTANT)) DEALLOCATE(HOST_PARAMS_CONSTANT) 
        if (ALLOCATED(HOST_PARAMS_CURRENT)) DEALLOCATE(HOST_PARAMS_CURRENT) 
        ALLOCATE(HOST_PARAMS_CURRENT(nparams))
        ALLOCATE(HOST_PARAMS_CONSTANT(nparams))
        HOST_PARAMS_CONSTANT = params
        HOST_PARAMS_CURRENT = params
        HOST_MODEL_NAME = TRIM(model_name)        
        HOST_STRUCTURE_SET = .TRUE.
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
        IF (.NOT. ALLOCATED(HOST_TIMES)) RETURN

        CALL query_current_kinematics(t)

    END SUBROUTINE update_hostcluster_state

    SUBROUTINE eval_force(nparticles, x, y, z, ax, ay, az)
        INTEGER, INTENT(IN) :: nparticles
        REAL*8, INTENT(IN), DIMENSION(nparticles) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(nparticles) :: ax, ay, az
        REAL*8, DIMENSION(nparticles,3) :: force_temp

        REAL*8, DIMENSION(nparticles) :: dx,dy,dz

        dx = x - HOST_X_CURRENT
        dy = y - HOST_Y_CURRENT
        dz = z - HOST_Z_CURRENT

        call model_force(nparticles,dx,dy,dz,force_temp)
        ax = force_temp(:,1)
        ay = force_temp(:,2)
        az = force_temp(:,3)

    END SUBROUTINE eval_force


    !! TO INTERFACE WITH SIMULATOR
    ! in hostcluster.f90 (inside CONTAINS)
    SUBROUTINE get_kinematics(ntimes, t, x, y, z, vx, vy, vz)
        INTEGER, INTENT(IN) :: ntimes
        REAL*8, INTENT(OUT), DIMENSION(ntimes) :: t, x, y, z, vx, vy, vz


        IF (.NOT. HOST_REGISTERED) THEN
            PRINT*, "WARNING: get_kinematics: host is not registered"
            RETURN
        END IF

        IF (.NOT. ALLOCATED(HOST_TIMES)) THEN
            PRINT*, "WARNING: get_kinematics: host kinematics are not set"
            RETURN
        END IF

        IF (SIZE(HOST_TIMES) /= ntimes) THEN
            PRINT*, "WARNING: get_kinematics: ntimes mismatch"
            RETURN
        END IF

        t  = HOST_TIMES
        x  = HOST_X
        y  = HOST_Y
        z  = HOST_Z
        vx = HOST_VX
        vy = HOST_VY
        vz = HOST_VZ

    END SUBROUTINE GET_KINEMATICS    

    ! extract the structural params
    subroutine get_structure(n_params,model_name,constant_params)
        INTEGER, INTENT(IN)                         :: n_params
        CHARACTER(LEN=64), INTENT(OUT)              :: model_name 
        REAL*8, INTENT(OUT), DIMENSION(n_params)    :: constant_params

        IF (.NOT. HOST_STRUCTURE_SET) THEN
            PRINT*, "WARNING: get_structure: host structure is not HOST_STRUCTURE_SET"
            RETURN
        END IF

        IF (HOST_NPARAMS /= n_params) THEN
            PRINT*, "WARNING: get_structure: n_params /= HOST_NPARAMS"
            RETURN
        END IF
        
        constant_params = HOST_PARAMS_CONSTANT
        model_name = HOST_MODEL_NAME

    END SUBROUTINE get_structure

    !!! MODELS 

    ! ANALYTICAL MODELS
    SUBROUTINE plummer_force(n, x, y, z, force)
        
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n,3) :: force
        REAL*8, DIMENSION(n) :: r, amod
        REAL*8 :: m, b

        m = HOST_PARAMS_CURRENT(1)
        b = HOST_PARAMS_CURRENT(2)
        r = SQRT(x*x + y*y + z*z)
        amod = -G_HOSTCLUSTER*m / (r*r + b*b)**1.5

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

        m = HOST_PARAMS_CURRENT(1)
        b = HOST_PARAMS_CURRENT(2)
        r = SQRT(x*x + y*y + z*z)
        phi = -G_HOSTCLUSTER*m / SQRT(r*r + b*b)
    END SUBROUTINE plummer_potential

END MODULE hostcluster
