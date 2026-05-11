MODULE simulator
    ! UX: 
    ! (1) set necessary values and physics module (order independent): 
    !   1a. simulator.set_initial_conditions(x,y,z,vx,vy,vz)
    !   1b. simulator.set_scheme("name", integration_params)
    !   1b. simulator.addcomponent("model", gravity_params)
    ! (2) add optional physics - (order independent)
    !   ex. host
    !       2bi.  simulator.init_host(t_host,x_host,v_host)
    !       2bii. simulator.set_host_model("plummer", params)
    ! (3) set i/o info - (order independent, can be set before anything else)
    !       simulator.write_snapshot(nskip,dir,fnamebase) - snapshot of all particles in one file
    !       simulator.write_orbits(nskip,dir,fnamebase) - each particle gets its own file and is continually appended
    ! (4) setbackward orbit (optional)
    !     simulator.setbackwardorbit()
    !       an optional call that will change the sign of the velocities and ensure the timestamps are descending  
    ! (5) finalize
    !       simulator.finalize()
    !           this subroutine performs all the proper checks before running the simulator
    ! (6) run()
    !       simulator.run()
    !           launches the computation based on all the entered information 
    ! (7) access the data
    !       simulator.x     ! The final positions
    !       simulator.xt    ! Some saved orbital trajectories based on the memory constraint
    ! (8) CLEAR()
    !       simulator.clear() the module before using again 
    IMPLICIT NONE 



    ! A derived type for orchestating the module 
    TYPE :: state_t
        ! REQUIRED
        LOGICAL :: initial_conditions_set = .FALSE.  ! the initial conditions have been entered
        LOGICAL :: gravity_finalized = .FALSE.       ! gravity module is configured and finalized
        LOGICAL :: scheme_set = .FALSE.              ! the integration technique is set 
        ! OPTIONAL PHYSICS
        LOGICAL :: host_enabled = .FALSE.                   ! OPTIONAL : using a host perturber 
        LOGICAL :: bar_enabled = .FALSE.                    ! OPTIONAL : using a galactic bar 
        LOGICAL :: perturbers_enabled = .FALSE.             ! OPTIONAL : using host perturbers
        ! TIME/SCHEME
        LOGICAL :: timestamps_from_user = .FALSE.
        LOGICAL :: backward_orbit = .FALSE.
        ! I/O 
        ! other state stuff
        LOGICAL :: finalized = .FALSE.
    END TYPE 

    ! Define scheme interface
    ABSTRACT INTERFACE
        SUBROUTINE scheme_step_interface()
        END SUBROUTINE scheme_step_interface
    END INTERFACE    

    TYPE(state_t), PRIVATE :: state

    ! numerical scheme varaiables
    CHARACTER(LEN=64), PRIVATE :: scheme_name = ""
    REAL*8, DIMENSION(:), ALLOCATABLE, PRIVATE :: scheme_params
    PROCEDURE(scheme_step_interface), POINTER, PRIVATE :: scheme => NULL()
    
    INTEGER, PUBLIC :: Nparticles
    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: x,y,z,vx,vy,vz
    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: timestamps
    INTEGER, PRIVATE :: NSTEPS, current_step


    CONTAINS 
   SUBROUTINE CLEAR()
    ! erases everything so it may be used again 
        if (state%initial_conditions_set) then 
            DEALLOCATE(x,y,z,vx,vy,vz)
            state%initial_conditions_set = .FALSE.
        END IF 
    END SUBROUTINE CLEAR

    SUBROUTINE finalize()
        ! builds tables, creates temp/dirs and runs last checks before computing.
        state%finalized = .TRUE.
    END subroutine finalize
    
    SUBROUTINE setinitialconditions(N,xin,yin,zin,vxin,vyin,vzin)
        INTEGER, INTENT(in) :: N 
        REAL*8, DIMENSION(N), INTENT(IN) :: xin,yin,zin,vxin,vyin,vzin
        Nparticles = N
        allocate(x(Nparticles),y(Nparticles),z(Nparticles),vx(Nparticles),vy(Nparticles),vz(Nparticles))
        x = xin
        y = yin 
        z = zin 
        vx = vxin 
        vy = vyin
        vz = vzin 
        state%initial_conditions_set = .TRUE.
    END SUBROUTINE setinitialconditions

    ! integration parameters
    SUBROUTINE build_fixed_timestamps()
        INTEGER :: i
        REAL*8 :: t0, dtmag, sgn

        IF (.NOT. ALLOCATED(scheme_params)) THEN
            PRINT*, "ERROR: scheme_params not set"
            RETURN
        END IF

        IF (SIZE(scheme_params) < 3) THEN
            PRINT*, "ERROR: fixed-step schemes require params=[t0, dt, nsteps]"
            RETURN
        END IF

        t0 = scheme_params(1)
        dtmag = ABS(scheme_params(2))
        nsteps = INT(scheme_params(3))

        IF (nsteps < 1) THEN
            PRINT*, "ERROR: nsteps must be >= 1"
            RETURN
        END IF

        IF (ALLOCATED(timestamps)) DEALLOCATE(timestamps)
        ALLOCATE(timestamps(nsteps + 1))

        IF (state%backward_orbit) THEN
            sgn = -1.0D0
        ELSE
            sgn =  1.0D0
        END IF

        DO i = 1, nsteps + 1
            timestamps(i) = t0 + sgn * dtmag * DBLE(i - 1)
        END DO
    END SUBROUTINE build_fixed_timestamps    

    SUBROUTINE setschema(name, params, timestamps_optional)
        CHARACTER(LEN=*), INTENT(IN) :: name
        REAL*8, DIMENSION(:), INTENT(IN) :: params 
        REAL*8, DIMENSION(:), INTENT(IN), OPTIONAL :: timestamps_optional

        IF (allocated(scheme_params)) DEALLOCATE(scheme_params)
        ALLOCATE(scheme_params(SIZE(params)))
        scheme_params = params
        scheme_name = TRIM(name)

        SELECT CASE (TRIM(NAME))
            CASE ("leapfrog")
                scheme=>leapfrog
            CASE ("forest_ruth")
                scheme=> forest_ruth
            CASE DEFAULT
                PRINT*, "ERROR: unknown scheme:, ", TRIM(name)
                NULLIFY(scheme)
                state%scheme_set = .false.
                RETURN 
        END SELECT 

        IF (PRESENT(timestamps_optional)) THEN 
            IF (SIZE(timestamps_optional) <2 ) then 
                print*, "ERROR timestamps must have at least 2 entries"
                state%scheme_set = .FALSE.
                RETURN 
            END IF 
            
            IF (ALLOCATED(timestamps)) DEALLOCATE(timestamps)
            ALLOCATE(timestamps(SIZE(timestamps_optional)))
            timestamps = timestamps_optional
            nsteps = SIZE(timestamps) - 1
            state%timestamps_from_user = .TRUE.
        
        ELSE 
            state%timestamps_from_user = .FALSE.
            call build_fixed_timestamps()
        END IF 

        state%scheme_set = .TRUE.

    END SUBROUTINE setschema

    ! THE NUMERICAL SCHEMES
    SUBROUTINE leapfrog()
        REAL*8, DIMENSION(Nparticles) :: xtemp, ytemp, ztemp
        REAL*8, DIMENSION(Nparticles) :: fx, fy, fz
        REAL*8 :: dt_step

        IF (current_step > nsteps) RETURN

        dt_step = timestamps(current_step + 1) - timestamps(current_step)

        xtemp = x + 0.5D0 * dt_step * vx
        ytemp = y + 0.5D0 * dt_step * vy
        ztemp = z + 0.5D0 * dt_step * vz

        ! Placeholder until force registry is wired.
        fx = 0.0D0
        fy = 0.0D0
        fz = 0.0D0

        vx = vx + dt_step * fx
        vy = vy + dt_step * fy
        vz = vz + dt_step * fz

        x = xtemp + 0.5D0 * dt_step * vx
        y = ytemp + 0.5D0 * dt_step * vy
        z = ztemp + 0.5D0 * dt_step * vz

        current_step = current_step + 1
    END SUBROUTINE leapfrog

    SUBROUTINE forest_ruth()

        print*, "to complete"

    END SUBROUTINE forest_ruth

    ! HELPER FUNCTIONS 
        LOGICAL FUNCTION is_strictly_increasing(t)
        REAL*8, DIMENSION(:), INTENT(IN) :: t
        INTEGER :: i
        is_strictly_increasing = .TRUE.
        DO i = 1, SIZE(t)-1
            IF (t(i+1) <= t(i)) THEN
                is_strictly_increasing = .FALSE.
                RETURN
            END IF
        END DO
    END FUNCTION is_strictly_increasing
    
    LOGICAL FUNCTION is_strictly_decreasing(t)
        REAL*8, DIMENSION(:), INTENT(IN) :: t
        INTEGER :: i
        is_strictly_decreasing = .TRUE.
        DO i = 1, SIZE(t)-1
            IF (t(i+1) >= t(i)) THEN
                is_strictly_decreasing = .FALSE.
                RETURN
            END IF
        END DO
    END FUNCTION is_strictly_decreasing    


END MODULE simulator
