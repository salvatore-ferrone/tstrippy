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
        LOGICAL :: orbits_allocated = .FALSE.
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

    ! for saving some trajectories
    REAL*8, PUBLIC :: orbit_ram_limit_MB = 1024.0D0 ! the default
    INTEGER, PRIVATE :: n_particles_orbit = 1 
    INTEGER, PRIVATE :: nskip_orbit_timestamps = 1 
    INTEGER, PRIVATE :: nvars_orbits = 7 ! time and phase space 
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: orbits

    ! some other limits
    REAL*8, PUBLIC :: max_ram_MB = 1024.0D0


    CONTAINS 


    SUBROUTINE CLEAR()
        ! Positions / velocities
        IF (ALLOCATED(x))  DEALLOCATE(x)
        IF (ALLOCATED(y))  DEALLOCATE(y)
        IF (ALLOCATED(z))  DEALLOCATE(z)
        IF (ALLOCATED(vx)) DEALLOCATE(vx)
        IF (ALLOCATED(vy)) DEALLOCATE(vy)
        IF (ALLOCATED(vz)) DEALLOCATE(vz)

        ! Timestamps and scheme params
        IF (ALLOCATED(timestamps))   DEALLOCATE(timestamps)
        IF (ALLOCATED(scheme_params)) DEALLOCATE(scheme_params)

        ! Scheme pointer and name
        NULLIFY(scheme)
        scheme_name = ""

        ! Counters
        Nparticles   = 0
        nsteps       = 0
        current_step = 1

        ! Reset all state flags to defaults
        state = state_t()

    END SUBROUTINE CLEAR


    SUBROUTINE compute_memory_particle_limit(particle_limit)
        INTEGER, INTENT(OUT) :: particle_limit
        INTEGER, PARAMETER :: nvariables = 6
        REAL*8, PARAMETER :: datasize_MB = 8.0D-6

        IF (max_ram_MB <= 0.0D0) THEN
            particle_limit = 0
            RETURN
        END IF

        particle_limit = INT(max_ram_MB / (DBLE(nvariables) * datasize_MB), kind=kind(particle_limit))
    END SUBROUTINE compute_memory_particle_limit

    SUBROUTINE setinitialconditions(N,xin,yin,zin,vxin,vyin,vzin)
        INTEGER, INTENT(in) :: N 
        INTEGER :: max_particles
        REAL*8, DIMENSION(N), INTENT(IN) :: xin,yin,zin,vxin,vyin,vzin
        ! Check if this violates configured particle-memory limit.

        IF (N < 1) THEN
            PRINT*, "ERROR in setinitialconditions: N must be >= 1"
            RETURN
        END IF

        CALL compute_memory_particle_limit(max_particles)
        if (N.gt.max_particles) then 
            PRINT*, "WARNING: Estimated memory for particles + overhead exceeds limit"
            print*, "   Max particles is", max_particles, "user entered: ", N
            print*, "   Break up the simulation or increase simulator.max_ram_MB"
            return 
        END IF 

        Nparticles = N
        IF (ALLOCATED(x))  DEALLOCATE(x)
        IF (ALLOCATED(y))  DEALLOCATE(y)
        IF (ALLOCATED(z))  DEALLOCATE(z)
        IF (ALLOCATED(vx)) DEALLOCATE(vx)
        IF (ALLOCATED(vy)) DEALLOCATE(vy)
        IF (ALLOCATED(vz)) DEALLOCATE(vz)
        allocate(x(Nparticles),y(Nparticles),z(Nparticles),vx(Nparticles),vy(Nparticles),vz(Nparticles))
        x = xin
        y = yin 
        z = zin 
        vx = vxin 
        vy = vyin
        vz = vzin 
        state%initial_conditions_set = .TRUE.
        state%finalized = .FALSE.
    END SUBROUTINE setinitialconditions

    SUBROUTINE setscheme(name, params, nparams)
        CHARACTER(LEN=*), INTENT(IN) :: name
        INTEGER, INTENT(IN) :: nparams
        REAL*8, DIMENSION(nparams), INTENT(IN) :: params 

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

        state%scheme_set = .TRUE.

    END SUBROUTINE setscheme

    SUBROUTINE settimestamps(tstamps, nt)
        INTEGER, INTENT(IN) :: nt
        REAL*8, DIMENSION(nt), INTENT(IN) :: tstamps
        IF (ALLOCATED(timestamps)) DEALLOCATE(timestamps)
        ALLOCATE(timestamps(nt))
        timestamps = tstamps
        nsteps = nt - 1
        state%timestamps_from_user = .TRUE.
        state%finalized = .FALSE.
    END SUBROUTINE settimestamps

    SUBROUTINE setbackwardorbit()
        state%backward_orbit = .TRUE.
        state%finalized = .FALSE.
    END SUBROUTINE setbackwardorbit

    SUBROUTINE run()
        
        call finalize()

        if (.not. state%finalized) then 
            print*, "ERROR: finalize failed. Cannot Run"
            RETURN 
        END IF 


    END SUBROUTINE run 

    SUBROUTINE finalize()
        LOGICAL :: should_return = .FALSE.
        state%finalized = .FALSE.

        IF (.NOT. state%initial_conditions_set) THEN
            PRINT*, "ERROR: initial conditions are not set"
            RETURN
        END IF

        IF (.NOT. ASSOCIATED(scheme)) THEN
            PRINT*, "ERROR: integration scheme not configured"
            should_return = .TRUE.
        END IF

        ! Resolve timestamps: build if user did not provide them
        IF (.NOT. state%timestamps_from_user) THEN
            CALL build_fixed_timestamps()
        END IF

        ! Universal checks on timestamps regardless of source
        IF (.NOT. ALLOCATED(timestamps)) THEN
            PRINT*, "ERROR: timestamps not available"
            should_return = .TRUE.
        END IF

        IF (SIZE(timestamps) < 2) THEN
            PRINT*, "ERROR: timestamps must have at least 2 entries"
            should_return = .TRUE.
        END IF

        IF (state%backward_orbit) THEN
            IF (.NOT. is_strictly_decreasing(timestamps)) THEN
                PRINT*, "ERROR: backward orbit requires strictly decreasing timestamps"
                should_return = .TRUE.
            END IF
        ELSE
            IF (.NOT. is_strictly_increasing(timestamps)) THEN
                PRINT*, "ERROR: forward orbit requires strictly increasing timestamps"
                should_return = .TRUE.
            END IF
        END IF

        if (should_return) return
        
        if (.not.state%orbits_allocated) THEN 
            print*, "ALLOCATING ORBITS!"
            CALL allocate_orbits(nskip_orbit_timestamps)
        END IF 
        

        nsteps = SIZE(timestamps) - 1
        current_step = 1
        state%scheme_set = .TRUE.
        state%finalized = .TRUE.

    END subroutine finalize
    
    ! CALLS MADE BY FINALIZE 
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

    subroutine trim_orbits(NSKIP)
        INTEGER, INTENT(IN) :: NSKIP 
        nskip_orbit_timestamps = NSKIP
    END SUBROUTINE trim_orbits

    SUBROUTINE allocate_orbits(NSKIP)
        INTEGER, INTENT(IN):: NSKIP
        LOGICAL :: QUIT = .False.
        REAL*8 :: memory_estimate
        REAL*8 :: memory_per_particle_per_step_MB = 8.0D-6
        INTEGER :: NSTEPS_ORBITS
        INTEGER :: NSAVED_ORBITS
        
        if (NSKIP.lt.1) then
            print*, "ERROR in allocate_orbits: NSKIP must be a positive integer"
            QUIT = .TRUE.
        END IF 

        IF (allocated(orbits)) DEALLOCATE(orbits)
        state%orbits_allocated = .FALSE.
        
        IF (.NOT.state%initial_conditions_set) then 
            print*, "ERROR in allocate_orbits: The initial conditions must be set before allocating space for the orbit trajectories"
            QUIT = .TRUE.
        END IF 

        IF (.NOT.state%scheme_set) then 
            print*, "ERROR in allocate_orbits: the scheme must be set before allocating space for the orbit trajectories"
            QUIT = .TRUE. 
        END IF 

        IF (Nparticles.lt.1) then
            print*, "ERROR in allocate_orbits: Nparticles must be >= 1"
            QUIT = .TRUE.
        END IF

        if (orbit_ram_limit_MB <= 0.0D0) THEN 
            print*, "ERROR in allocate_orbits: orbit_ram_limit_MB must be > 0"
            n_particles_orbit = 0
            QUIT = .TRUE.
        END IF 

        if (QUIT) THEN 
            RETURN 
        END IF 

        NSTEPS_ORBITS = NSTEPS / NSKIP
        NSAVED_ORBITS = NSTEPS_ORBITS + 1
        
        ! Estimate storage for all particles over all saved orbit snapshots.
        memory_estimate = DBLE(Nparticles) * DBLE(NSAVED_ORBITS) * DBLE(nvars_orbits) * memory_per_particle_per_step_MB

        IF (memory_estimate .gt. orbit_ram_limit_MB) THEN
            n_particles_orbit = INT(orbit_ram_limit_MB  / (DBLE(nvars_orbits) * DBLE(NSAVED_ORBITS) * memory_per_particle_per_step_MB), &
                kind=kind(n_particles_orbit) )
        ELSE
            n_particles_orbit = Nparticles
        END IF

        n_particles_orbit = max(0, min(n_particles_orbit, Nparticles))

        ! check if no particles will be saved 
        if (n_particles_orbit.lt.1) THEN 
            print*, "WARNING. No orbital trajectories will be allocated."
            print*, "   Memory allocation limit:", orbit_ram_limit_MB, "MBytes"
            print*, "   Size of full simulation:", memory_estimate, "MBytes"
            print*, "   Size of a single orbit:", int(memory_estimate/DBLE(Nparticles))
            print*, "   To get a sample of orbits, call `simulator.trim_orbits(NSKIP)`, to skip over intermediate timesteps"
        ELSE IF (memory_estimate .gt. orbit_ram_limit_MB) THEN 
            print*, "NOTE TO USER."
            print*, "   Only subset of the orbital trajectories will be allocated to `simulator.orbits`"
            print*, "   ALLOCATING: ", n_particles_orbit, "/", Nparticles, "particles"
            print*, "   IF all orbits are desired, reduce total integration time or call `trim_orbits(NSKIP)` to skip over intermediate timesteps"
        end IF 

        allocate(orbits(NSAVED_ORBITS,nvars_orbits,n_particles_orbit))
        state%orbits_allocated=.TRUE.

    END SUBROUTINE allocate_orbits

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
