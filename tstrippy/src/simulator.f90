MODULE simulator
    USE gravity, ONLY: gravity_clear => clear, &
                       gravity_set_gravitational_constant => set_gravitational_constant, &
                       gravity_add_component => add_component, &
                       gravity_add_component_agama => add_component_agama, &
                       gravity_add_component_agama_from_file => add_component_agama_from_file, &
                       gravity_finalize => finalize, &
                       gravity_force => force, &
                       gravity_potential => potential, &
                       GRAVITY_FINALIZED, &
                       GRAVITY_NCOMP
    USE hostcluster, ONLY: hostcluster_clear                => clear, &
                           hostcluster_add                  => add_hostcluster, &
                           hostcluster_configure_kinematics => configure_hostcluster_kinematics, &
                           hostcluster_configure_model      => configure_hostcluster_structure, &
                           hostcluster_finalize             => finalize_hostcluster, &
                           hostcluster_update_state         => update_hostcluster_state, &
                           hostcluster_eval_force           => eval_force, &
                           hostcluster_set_gravitational_constant => set_gravitational_constant,&
                           hostcluster_get_kinematics       => get_kinematics,&
                           hostcluster_get_structure        => get_structure,&
                           hostcluster_configure_structure_parameter_table  => configure_hostcluster_structure_parameter_table, &
                           hostcluster_initialize_ionization_state          => initialize_ionization_state, &
                           hostcluster_update_ionization_state              => update_ionization_state, &
                           hostcluster_get_ionization_state                 => get_ionization_state, &
                           HOST_REGISTERED, &
                           HOST_FINALIZED
    use flyingspheres, ONLY: flyingspheres_clear => clear,&
                                flyingspheres_initialize => initialize_flyingspheres,&
                                flyingsphere_set => set_flyingsphere,&
                                flyingspheres_finalize => finalize_flyingspheres,&
                                flyingspheres_update_state => update_state,&
                                flyingspheres_force => eval_total_force,&
                                FLYINGSPHERES_REGISTERED,&
                                FLYINGSPHERES_FINALIZED

    USE mathutils, ONLY: is_strictly_increasing, is_strictly_decreasing
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
    TYPE, PRIVATE :: state_t
        ! REQUIRED
        LOGICAL :: initial_conditions_set = .FALSE.  ! the initial conditions have been entered
        LOGICAL :: gravity_finalized = .FALSE.       ! gravity module is configured and finalized
        LOGICAL :: scheme_set = .FALSE.              ! the integration technique is set 
        ! OPTIONAL PHYSICS
        LOGICAL :: host_enabled = .FALSE.                   ! OPTIONAL : using a host perturber 
        LOGICAL :: bar_enabled = .FALSE.                    ! OPTIONAL : using a galactic bar 
        LOGICAL :: flyingspheres_enabled = .FALSE.             ! OPTIONAL : using host perturbers
        ! TIME/SCHEME
        LOGICAL :: timestamps_from_user = .FALSE.
        LOGICAL :: backward_orbit = .FALSE.
        ! I/O 
        LOGICAL :: orbits_allocated = .FALSE.
        LOGICAL :: writesnapshots = .FALSE.
        LOGICAL :: writeorbits = .FALSE.
        ! other state stuff
        LOGICAL :: finalized = .FALSE.
        LOGICAL :: DONTRUN = .FALSE.
        LOGICAL :: run_success = .FALSE.
    END TYPE 

    ! Define scheme interface
    ABSTRACT INTERFACE
        SUBROUTINE scheme_step_interface()
        END SUBROUTINE scheme_step_interface

        SUBROUTINE force_provider_interface(t, n, x, y, z, ax, ay, az)
            REAL*8, INTENT(IN) :: t
            INTEGER, INTENT(IN) :: n
            REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
            REAL*8, INTENT(OUT), DIMENSION(n) :: ax, ay, az
        END SUBROUTINE force_provider_interface
    END INTERFACE    

    TYPE, PRIVATE :: force_provider_t
        CHARACTER(LEN=64) :: name = ""
        PROCEDURE(force_provider_interface), POINTER, NOPASS :: eval => NULL()
    END TYPE force_provider_t

    TYPE(state_t), PRIVATE :: state

    LOGICAL :: FUN = .False.
    ! numerical scheme varaiables
    CHARACTER(LEN=64), PRIVATE :: scheme_name = ""
    REAL*8, DIMENSION(:), ALLOCATABLE, PRIVATE :: scheme_params
    PROCEDURE(scheme_step_interface), POINTER, PRIVATE :: scheme => NULL()
    
    INTEGER, PUBLIC :: Nparticles
    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: x,y,z,vx,vy,vz
    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: x_initial, y_initial, z_initial
    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: vx_initial, vy_initial, vz_initial
    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: timestamps
    INTEGER, PRIVATE :: NSTEPS, current_step
    REAL*8, PUBLIC :: currenttime
    
    REAL*8, PRIVATE :: yoshida_w, yoshida_c1, yoshida_c2, yoshida_c3, yoshida_c4
    REAL*8, PRIVATE :: yoshida_d1, yoshida_d2, yoshida_d3, yoshida_d4    

    INTEGER, PARAMETER, PRIVATE :: MAX_ACTIVE_FORCE_PROVIDERS = 8
    TYPE(force_provider_t), DIMENSION(MAX_ACTIVE_FORCE_PROVIDERS), PRIVATE :: active_force_providers
    INTEGER, PRIVATE :: nactive_force_providers = 0

    ! for saving some trajectories
    REAL*8, PARAMETER :: DEFAULT_ORBIT_RAM_LIMIT_MB = 1024.0D0
    REAL*8, PARAMETER :: DEFAULT_MAX_RAM_MB = 1024.0D0

    REAL*8, PUBLIC :: orbit_ram_limit_MB = DEFAULT_ORBIT_RAM_LIMIT_MB
    INTEGER, PRIVATE :: n_particles_orbit = 1 
    INTEGER, PRIVATE :: nskip_orbit_timestamps = 1 
    INTEGER, PRIVATE :: nvars_orbits = 6 ! phase space 
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: orbits
    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: orbits_timestamps
    
    ! persistent file handles for direct orbit writing (no per-step open/close)
    INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: orbit_file_units

    INTEGER, PARAMETER :: DEFAULT_FILEUNITBASE_WRITESNAPSHOTS = 12345
    INTEGER, PUBLIC :: FILEUNITBASE_WRITESNAPSHOTS = DEFAULT_FILEUNITBASE_WRITESNAPSHOTS
    INTEGER, PUBLIC :: FILEUNITBASE_WRITEORBITS    = DEFAULT_FILEUNITBASE_WRITESNAPSHOTS + 1
    INTEGER, PUBLIC :: nskip_writesnapshots = 1
    INTEGER, PUBLIC :: nskip_writeorbits    = 1
    CHARACTER(LEN=500), PUBLIC :: directory_snapshots = ""
    CHARACTER(LEN=500), PUBLIC :: basename_snapshots  = ""
    CHARACTER(LEN=500), PUBLIC :: directory_orbits    = ""
    CHARACTER(LEN=500), PUBLIC :: basename_orbits     = ""
    ! some other limits
    REAL*8, PUBLIC :: max_ram_MB = DEFAULT_MAX_RAM_MB
    ! flat public state for Python inspection after run()
    LOGICAL, PUBLIC :: RUN_SUCCESS         = .FALSE.
    LOGICAL, PUBLIC :: DID_WRITE_SNAPSHOTS = .FALSE.
    LOGICAL, PUBLIC :: DID_WRITE_ORBITS    = .FALSE.
    LOGICAL, PUBLIC :: BACKWARD_ORBIT_ENABLED = .FALSE.
    CHARACTER(LEN=64), PUBLIC :: SCHEME_METHOD = ""
    INTEGER, PUBLIC :: SCHEME_NPARAMS = 0
    REAL*8, DIMENSION(16), PUBLIC :: SCHEME_PARAMETERS = 0.0D0
    ! Timing outputs (all in seconds)
    REAL*8, PUBLIC :: TIMER_RUN_SECONDS = 0.0D0
    REAL*8, PUBLIC :: TIMER_FINALIZE_SECONDS = 0.0D0
    REAL*8, PUBLIC :: TIMER_SCHEME_SECONDS = 0.0D0
    REAL*8, PUBLIC :: TIMER_WRITE_SNAPSHOTS_SECONDS = 0.0D0
    REAL*8, PUBLIC :: TIMER_WRITE_ORBITS_SECONDS = 0.0D0
    REAL*8, PUBLIC :: TIMER_FLYINGSPHERES_SECONDS = 0.0D0
    ! Private timing state (for run)
    INTEGER(KIND=8), PRIVATE :: c_run_start, c_run_end, c_scheme_accum, c_scheme_start
    INTEGER(KIND=8), PRIVATE :: c_write_snap_accum, c_write_orb_accum
    INTEGER(KIND=8), PRIVATE :: c_flyingspheres_start = 0 
    INTEGER(KIND=8), PRIVATE :: c_flyingspheres_end = 0 
    INTEGER(KIND=8), PRIVATE :: c_flyingspheres_accum = 0 
    INTEGER(KIND=8), PRIVATE :: rate_clock, cmax_clock

    ! PUBLIC :: simulator_cleargravitycomponents, simulator_set_gravitational_constant
    ! PUBLIC :: simulator_add_component, simulator_finalizegravity
    ! PUBLIC :: simulator_force, simulator_potential

    INTEGER, PUBLIC :: N_HOST_ORBIT_TIME_STAMPS = 0
    INTEGER, PUBLIC :: N_HOST_STRUCTURE_PARAMETERS = 0 


    CONTAINS 
    
    !!!!! CALLS WHERE THE USER INTERFACES WITH THE MODULE
    SUBROUTINE set_gravitational_constant(g)
        REAL*8, INTENT(IN) :: g

        CALL hostcluster_set_gravitational_constant(g)
        CALL gravity_set_gravitational_constant(g)
        CALL clear_force_registry()
        state%gravity_finalized = GRAVITY_FINALIZED
        state%finalized = .FALSE.
    END SUBROUTINE set_gravitational_constant

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

        IF (ALLOCATED(x_initial))  DEALLOCATE(x_initial)
        IF (ALLOCATED(y_initial))  DEALLOCATE(y_initial)
        IF (ALLOCATED(z_initial))  DEALLOCATE(z_initial)
        IF (ALLOCATED(vx_initial)) DEALLOCATE(vx_initial)
        IF (ALLOCATED(vy_initial)) DEALLOCATE(vy_initial)
        IF (ALLOCATED(vz_initial)) DEALLOCATE(vz_initial)
        ALLOCATE(x_initial(Nparticles), y_initial(Nparticles), z_initial(Nparticles))
        ALLOCATE(vx_initial(Nparticles), vy_initial(Nparticles), vz_initial(Nparticles))
        x_initial = xin
        y_initial = yin
        z_initial = zin
        vx_initial = vxin
        vy_initial = vyin
        vz_initial = vzin

        IF (ALLOCATED(orbits)) DEALLOCATE(orbits)
        IF (ALLOCATED(orbits_timestamps)) DEALLOCATE(orbits_timestamps)
        state%orbits_allocated = .FALSE.
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
        SCHEME_METHOD = TRIM(name)
        SCHEME_PARAMETERS = 0.0D0
        SCHEME_NPARAMS = MIN(nparams, SIZE(SCHEME_PARAMETERS))
        IF (SCHEME_NPARAMS > 0) THEN
            SCHEME_PARAMETERS(1:SCHEME_NPARAMS) = params(1:SCHEME_NPARAMS)
        END IF

        SELECT CASE (TRIM(NAME))
            CASE ("leapfrog")
                scheme=>leapfrog
            CASE ("forest_ruth")
                CALL compute_yoshida_coefficients()
                scheme=> forest_ruth
            CASE DEFAULT
                PRINT*, "ERROR: unknown scheme:, ", TRIM(name)
                NULLIFY(scheme)
                state%scheme_set = .false.
                SCHEME_METHOD = ""
                SCHEME_NPARAMS = 0
                SCHEME_PARAMETERS = 0.0D0
                RETURN 
        END SELECT

        IF (ALLOCATED(orbits)) DEALLOCATE(orbits)
        IF (ALLOCATED(orbits_timestamps)) DEALLOCATE(orbits_timestamps)
        state%orbits_allocated = .FALSE.
        state%scheme_set = .TRUE.
        state%finalized = .FALSE.

    END SUBROUTINE setscheme

    SUBROUTINE settimestamps(tstamps, nt)
        INTEGER, INTENT(IN) :: nt
        REAL*8, DIMENSION(nt), INTENT(IN) :: tstamps
        IF (ALLOCATED(timestamps)) DEALLOCATE(timestamps)
        ALLOCATE(timestamps(nt))
        timestamps = tstamps
        nsteps = nt - 1
        IF (ALLOCATED(orbits)) DEALLOCATE(orbits)
        IF (ALLOCATED(orbits_timestamps)) DEALLOCATE(orbits_timestamps)
        state%orbits_allocated = .FALSE.
        state%timestamps_from_user = .TRUE.
        state%finalized = .FALSE.
    END SUBROUTINE settimestamps

    SUBROUTINE setbackwardorbit()
        state%backward_orbit = .TRUE.
        BACKWARD_ORBIT_ENABLED = .TRUE.
        state%finalized = .FALSE.
    END SUBROUTINE setbackwardorbit

    subroutine trim_orbits(NSKIP)
        INTEGER, INTENT(IN) :: NSKIP 
        IF (NSKIP < 1) THEN
            PRINT*, "ERROR in trim_orbits: NSKIP must be >= 1"
            RETURN
        END IF

        nskip_orbit_timestamps = NSKIP
        IF (ALLOCATED(orbits)) DEALLOCATE(orbits)
        IF (ALLOCATED(orbits_timestamps)) DEALLOCATE(orbits_timestamps)
        state%orbits_allocated = .FALSE.
        state%finalized = .FALSE.
    END SUBROUTINE trim_orbits
    
    SUBROUTINE clear()
        ! Positions / velocities
        IF (ALLOCATED(x))  DEALLOCATE(x)
        IF (ALLOCATED(y))  DEALLOCATE(y)
        IF (ALLOCATED(z))  DEALLOCATE(z)
        IF (ALLOCATED(vx)) DEALLOCATE(vx)
        IF (ALLOCATED(vy)) DEALLOCATE(vy)
        IF (ALLOCATED(vz)) DEALLOCATE(vz)
        IF (ALLOCATED(x_initial))  DEALLOCATE(x_initial)
        IF (ALLOCATED(y_initial))  DEALLOCATE(y_initial)
        IF (ALLOCATED(z_initial))  DEALLOCATE(z_initial)
        IF (ALLOCATED(vx_initial)) DEALLOCATE(vx_initial)
        IF (ALLOCATED(vy_initial)) DEALLOCATE(vy_initial)
        IF (ALLOCATED(vz_initial)) DEALLOCATE(vz_initial)

        ! Timestamps and scheme params
        IF (ALLOCATED(timestamps))   DEALLOCATE(timestamps)
        IF (ALLOCATED(scheme_params)) DEALLOCATE(scheme_params)
        IF (ALLOCATED(orbits)) DEALLOCATE(orbits)
        IF (ALLOCATED(orbits_timestamps)) DEALLOCATE(orbits_timestamps)

        ! Scheme pointer and name
        NULLIFY(scheme)
        scheme_name = ""

        ! Counters
        Nparticles   = 0
        nsteps       = 0
        current_step = 1
        currenttime  = 0.0D0
        n_particles_orbit = 1
        nskip_orbit_timestamps = 1
        orbit_ram_limit_MB = DEFAULT_ORBIT_RAM_LIMIT_MB
        max_ram_MB = DEFAULT_MAX_RAM_MB

        ! I/O settings
        nskip_writesnapshots = 1
        nskip_writeorbits = 1
        directory_snapshots = ""
        basename_snapshots = ""
        directory_orbits = ""
        basename_orbits = ""
        
        ! Close and deallocate orbit file units
        IF (ALLOCATED(orbit_file_units)) THEN
            DO n_particles_orbit = 1, SIZE(orbit_file_units)
                IF (orbit_file_units(n_particles_orbit) > 0) THEN
                    CLOSE(orbit_file_units(n_particles_orbit))
                END IF
            END DO
            DEALLOCATE(orbit_file_units)
        END IF

        ! Reset flat public flags
        RUN_SUCCESS         = .FALSE.
        DID_WRITE_SNAPSHOTS = .FALSE.
        DID_WRITE_ORBITS    = .FALSE.
        BACKWARD_ORBIT_ENABLED = .FALSE.
        SCHEME_METHOD = ""
        SCHEME_NPARAMS = 0
        SCHEME_PARAMETERS = 0.0D0
        ! Reset timers
        TIMER_RUN_SECONDS = 0.0D0
        TIMER_FINALIZE_SECONDS = 0.0D0
        TIMER_SCHEME_SECONDS = 0.0D0
        TIMER_WRITE_SNAPSHOTS_SECONDS = 0.0D0
        TIMER_WRITE_ORBITS_SECONDS = 0.0D0
        TIMER_FLYINGSPHERES_SECONDS = 0.0D0
        c_run_start = 0
        c_run_end = 0
        c_scheme_accum = 0
        c_write_snap_accum = 0
        c_write_orb_accum = 0
        c_flyingspheres_accum = 0 
        c_flyingspheres_start = 0
        c_flyingspheres_end = 0
        ! Reset all state flags to defaults
        state = state_t()
        CALL clear_force_registry()
        ! clear the modules
        call hostcluster_clear()
        call gravity_clear()
        call flyingspheres_clear()
        N_HOST_ORBIT_TIME_STAMPS = 0 

    END SUBROUTINE clear

    SUBROUTINE finalize()
        LOGICAL :: should_return
        INTEGER(KIND=8) :: c_fin_start, c_fin_end
        should_return = .FALSE.
        state%finalized = .FALSE.

        CALL system_clock(c_fin_start, rate_clock, cmax_clock)

        IF (.NOT. state%initial_conditions_set) THEN
            PRINT*, "ERROR: initial conditions are not set"
            RETURN
        END IF

        IF (.NOT. ASSOCIATED(scheme)) THEN
            PRINT*, "ERROR: integration scheme not configured"
            should_return = .TRUE.
        END IF

        IF (GRAVITY_NCOMP < 1) THEN
            PRINT*, "ERROR: gravity is not configured. GRAVITY_NCOMP", GRAVITY_NCOMP
            should_return = .TRUE.
        ELSE IF (.NOT. GRAVITY_FINALIZED) THEN
            CALL gravity_finalize()
            IF (.NOT. GRAVITY_FINALIZED) THEN
                PRINT*, "ERROR: gravity finalize failed"
                should_return = .TRUE.
            END IF
        END IF

        state%gravity_finalized = GRAVITY_FINALIZED

        IF (state%host_enabled .OR. HOST_REGISTERED) THEN
            ! note. this should never trigger. However, it's important to note that 
            !   it doesn't make sense to use a host cluster without also using particles
            if (.not.state%initial_conditions_set) then 
                print*, "a host cluster can only be used when initial conditions are registered"
                should_return = .TRUE.
            END IF 

            CALL hostcluster_initialize_ionization_state(Nparticles)
            state%host_enabled = HOST_REGISTERED
            
            IF (.NOT. HOST_FINALIZED) THEN
                CALL hostcluster_finalize()
            END IF

            IF (.NOT. HOST_FINALIZED) THEN
                PRINT*, "ERROR: hostcluster finalize failed"
                should_return = .TRUE.
            END IF
            ! set the ionization state

        END IF

        IF (FLYINGSPHERES_REGISTERED) THEN
            state%flyingspheres_enabled = .TRUE.
            IF (.NOT. FLYINGSPHERES_FINALIZED) THEN
                CALL flyingspheres_finalize()
            END IF
            IF (.NOT. FLYINGSPHERES_FINALIZED) THEN
                PRINT*, "ERROR: flyingspheres finalize failed"
                should_return = .TRUE.
            END IF
        ELSE
            state%flyingspheres_enabled = .FALSE.
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

        nsteps = SIZE(timestamps) - 1
        current_step = 1
        currenttime = timestamps(current_step)

        if (should_return) return

        CALL rebuild_force_registry()
        IF (nactive_force_providers < 1) THEN
            PRINT*, "ERROR: no active force providers are registered"
            RETURN
        END IF

        if (.not.state%orbits_allocated) THEN
            CALL allocate_orbits(nskip_orbit_timestamps)
            IF (.NOT. state%orbits_allocated) THEN
                PRINT*, "ERROR: finalize failed because orbit allocation did not complete"
                RETURN
            END IF
        END IF

        state%scheme_set = .TRUE.
        state%finalized = .TRUE.

        CALL system_clock(c_fin_end, rate_clock, cmax_clock)
        IF (c_fin_end >= c_fin_start) THEN
            TIMER_FINALIZE_SECONDS = DBLE(c_fin_end - c_fin_start) / DBLE(rate_clock)
        ELSE
            TIMER_FINALIZE_SECONDS = DBLE((cmax_clock - c_fin_start) + c_fin_end + 1_8) / DBLE(rate_clock)
        END IF

    END subroutine finalize

    SUBROUTINE run()
        INTEGER :: istep, iorbit
        INTEGER(KIND=8) :: c_scheme_end, c_snap_start, c_snap_end, c_orb_start, c_orb_end

        CALL system_clock(c_run_start, rate_clock, cmax_clock)
        c_scheme_accum = 0
        c_write_snap_accum = 0
        c_write_orb_accum = 0
        c_flyingspheres_accum = 0 

        CALL finalize()

        IF (.NOT. state%finalized) THEN
            PRINT*, "ERROR: finalize failed. Cannot Run"
            RETURN
        END IF
        
        if (state%DONTRUN) then
            print*, "---------"
            print*, "ERROR IN simulator.RUN()"
            if (FUN) then
                print*, "  ________    _____      _____  ___________ "
                print*, " /  _____/   /  _  \    /     \ \_   _____/ "
                print*, "/   \  ___  /  /_\  \  /  \ /  \ |    __)_  "
                print*, "\    \_\  \/    |    \/    Y    \|        \ "
                print*, " \______  /\____|__  /\____|__  /_______  / "
                print*, "        \/         \/         \/        \/  "
                print*, "____________   _________________________    "
                print*, "\_____  \   \ /   /\_   _____/\______   \   "
                print*, " /   |   \   Y   /  |    __)_  |       _/   "
                print*, "/    |    \     /   |        \ |    |   \   "
                print*, "\_______  /\___/   /_______  / |____|_  /   "
                print*, "        \/                 \/         \/    "
            end if 
            print*, "---------"
            print*, "The guard DONTRUN is True. Aborting run!"
            return 
        end if 
        
        ! Open all orbit files before loop if writing orbits (one per particle, unlimited)
        IF (state%writeorbits .AND. Nparticles > 0) THEN
            CALL open_orbit_files()
            ! Write initial conditions immediately
            CALL write_orbit_records()
        END IF

        ! Save initial conditions as the first orbit snapshot
        iorbit = 1
        IF (n_particles_orbit > 0) THEN
            orbits_timestamps(iorbit) = timestamps(1)
            orbits(iorbit, 1, :) = x(1:n_particles_orbit)
            orbits(iorbit, 2, :) = y(1:n_particles_orbit)
            orbits(iorbit, 3, :) = z(1:n_particles_orbit)
            orbits(iorbit, 4, :) = vx(1:n_particles_orbit)
            orbits(iorbit, 5, :) = vy(1:n_particles_orbit)
            orbits(iorbit, 6, :) = vz(1:n_particles_orbit)
            iorbit = iorbit + 1
        END IF

        ! Main integration loop
        DO istep = 1, nsteps
            CALL system_clock(c_scheme_start)
            CALL scheme()  ! advances positions/velocities and increments current_step
            CALL system_clock(c_scheme_end)
            c_scheme_accum = c_scheme_accum + (c_scheme_end - c_scheme_start)

            ! ADD NOTE
            ! if there is a host cluster, check who is bound
            IF (state%host_enabled) then
                CALL hostcluster_update_state(currenttime)
                CALL hostcluster_update_ionization_state(Nparticles,x,y,z,vx,vy,vz)
            END IF 

            ! Save orbit snapshot every nskip_orbit_timestamps steps
            IF (n_particles_orbit > 0) THEN
                IF (MOD(istep, nskip_orbit_timestamps) == 0) THEN
                    orbits_timestamps(iorbit) = timestamps(current_step)
                    orbits(iorbit, 1, :) = x(1:n_particles_orbit)
                    orbits(iorbit, 2, :) = y(1:n_particles_orbit)
                    orbits(iorbit, 3, :) = z(1:n_particles_orbit)
                    orbits(iorbit, 4, :) = vx(1:n_particles_orbit)
                    orbits(iorbit, 5, :) = vy(1:n_particles_orbit)
                    orbits(iorbit, 6, :) = vz(1:n_particles_orbit)
                    iorbit = iorbit + 1
                END IF
            END IF

            ! Write snapshot if enabled
            IF (state%writesnapshots) THEN
                IF (MOD(istep, nskip_writesnapshots) == 0) THEN
                    CALL system_clock(c_snap_start)
                    CALL write_snapshot_file(istep)
                    CALL system_clock(c_snap_end)
                    c_write_snap_accum = c_write_snap_accum + (c_snap_end - c_snap_start)
                END IF
            END IF

            ! Write orbits directly (per-particle files, no memory limit)
            IF (state%writeorbits .AND. ALLOCATED(orbit_file_units)) THEN
                IF (MOD(istep, nskip_writeorbits) == 0) THEN
                    CALL system_clock(c_orb_start)
                    CALL write_orbit_records()
                    CALL system_clock(c_orb_end)
                    c_write_orb_accum = c_write_orb_accum + (c_orb_end - c_orb_start)
                END IF
            END IF
        END DO
        
        ! Close all orbit files at end
        IF (ALLOCATED(orbit_file_units)) THEN
            CALL close_orbit_files()
        END IF
        state%run_success = .TRUE.
        RUN_SUCCESS = .TRUE.
        DID_WRITE_SNAPSHOTS = state%writesnapshots
        DID_WRITE_ORBITS    = state%writeorbits
        if (FUN) then
            print*, ""
            print*, ""
            print*, "  /$$$$$$  /$$   /$$  /$$$$$$  /$$$$$$$$  /$$$$$$   /$$$$$$ "
            print*, " /$$__  $$| $$  | $$ /$$__  $$| $$_____/ /$$__  $$ /$$__  $$"
            print*, "| $$  \__/| $$  | $$| $$  \__/| $$      | $$  \__/| $$  \__/"
            print*, "|  $$$$$$ | $$  | $$| $$      | $$$$$   |  $$$$$$ |  $$$$$$ "
            print*, " \____  $$| $$  | $$| $$      | $$__/    \____  $$ \____  $$"
            print*, " /$$  \ $$| $$  | $$| $$    $$| $$       /$$  \ $$ /$$  \ $$"
            print*, "|  $$$$$$/|  $$$$$$/|  $$$$$$/| $$$$$$$$|  $$$$$$/|  $$$$$$/"
            print*, " \______/  \______/  \______/ |________/ \______/  \______/ "
            print*, ""
        end if 
        ! print*, "_________________________________________.____________________________.___."
        ! print*, "\__    ___/   _____/\__    ___/\______   \   \______   \______   \__  |   |"
        ! print*, "  |    |  \_____  \   |    |    |       _/   ||     ___/|     ___//   |   |"
        ! print*, "  |    |  /        \  |    |    |    |   \   ||    |    |    |    \____   |"
        ! print*, "  |____| /_______  /  |____|    |____|_  /___||____|    |____|    / ______|"
        ! print*, "                 \/                    \/                         \/       "
        ! print*, "___________________ _________________________   _________________________  "
        ! print*, "\_   _____/\_____  \\______   \_   _____/\   \ /   /\_   _____/\______   \ "
        ! print*, " |    __)   /   |   \|       _/|    __)_  \   Y   /  |    __)_  |       _/ "
        ! print*, " |     \   /    |    \    |   \|        \  \     /   |        \ |    |   \ "
        ! print*, " \___  /   \_______  /____|_  /_______  /   \___/   /_______  / |____|_  / "
        ! print*, "     \/            \/       \/        \/                    \/         \/  "

        ! Finalize timing measurements
        CALL system_clock(c_run_end)

        ! Calculate elapsed times (handle clock wraparound)
        IF (c_run_end >= c_run_start) THEN
            TIMER_RUN_SECONDS = DBLE(c_run_end - c_run_start) / DBLE(rate_clock)
        ELSE
            TIMER_RUN_SECONDS = DBLE((cmax_clock - c_run_start) + c_run_end + 1_8) / DBLE(rate_clock)
        END IF

        TIMER_SCHEME_SECONDS = DBLE(c_scheme_accum) / DBLE(rate_clock)
        TIMER_WRITE_SNAPSHOTS_SECONDS = DBLE(c_write_snap_accum) / DBLE(rate_clock)
        TIMER_WRITE_ORBITS_SECONDS = DBLE(c_write_orb_accum) / DBLE(rate_clock)
        TIMER_FLYINGSPHERES_SECONDS = DBLE(c_flyingspheres_accum) / DBLE(rate_clock)
        
    END SUBROUTINE run

    !!! INTERACTING WITH SPECIFIC MODULES
    
    !!! THE FLYING SPHERES 
    SUBROUTINE initialize_flyingspheres(n)
        integer, intent(in) :: n
        call flyingspheres_initialize(n)
        CALL clear_force_registry()
        state%flyingspheres_enabled = FLYINGSPHERES_REGISTERED
        state%finalized = .FALSE.
    END SUBROUTINE initialize_flyingspheres

    subroutine set_flyingsphere(i,model_name, nparams, structural_parameters, ntimestamps, txyz)
        INTEGER, INTENT(IN) :: i, nparams, ntimestamps
        character(len=64), INTENT(IN):: model_name
        REAl*8, DIMENSION(nparams), INTENT(IN) :: structural_parameters
        REAL*8, DIMENSION(4,ntimestamps),INTENT(IN) :: txyz
        call flyingsphere_set(i,model_name, nparams, structural_parameters, ntimestamps, txyz)
        CALL clear_force_registry()
        state%flyingspheres_enabled = FLYINGSPHERES_REGISTERED
        state%finalized = .FALSE.
    end subroutine set_flyingsphere

    SUBROUTINE finalize_flyingspheres()
        CALL flyingspheres_finalize()
        CALL clear_force_registry()
        state%flyingspheres_enabled = FLYINGSPHERES_REGISTERED .AND. FLYINGSPHERES_FINALIZED
        state%finalized = .FALSE.
    END SUBROUTINE finalize_flyingspheres
    
    !!! THE GRAVITY MODULE 
    SUBROUTINE cleargravitycomponents()
        CALL gravity_clear()
        CALL clear_force_registry()
        state%gravity_finalized = .FALSE.
        state%finalized = .FALSE.
    END SUBROUTINE cleargravitycomponents    

    SUBROUTINE add_component(model_name, params, nparams)
        CHARACTER(LEN=*), INTENT(IN) :: model_name
        INTEGER, INTENT(IN) :: nparams
        REAL*8, DIMENSION(nparams), INTENT(IN) :: params

        CALL gravity_add_component(model_name, params, nparams)
        CALL clear_force_registry()
        state%gravity_finalized = GRAVITY_FINALIZED
        state%finalized = .FALSE.
    END SUBROUTINE add_component

    SUBROUTINE add_component_agama(params)
        CHARACTER(LEN=*), INTENT(IN) :: params
        CALL gravity_add_component_agama(params)
        CALL clear_force_registry()
        state%gravity_finalized = GRAVITY_FINALIZED
        state%finalized = .FALSE.
    END SUBROUTINE add_component_agama   
    
    SUBROUTINE add_component_agama_from_file(inifilename)
        CHARACTER(LEN=*), INTENT(IN) :: inifilename
        CALL gravity_add_component_agama_from_file(inifilename)
        CALL clear_force_registry()
        state%gravity_finalized = GRAVITY_FINALIZED
        state%finalized = .FALSE.
    END SUBROUTINE add_component_agama_from_file        

    !!!!! THE HOST CLUSTER MODULE
    SUBROUTINE add_hostcluster()
        CALL hostcluster_add()
        CALL clear_force_registry()
        state%host_enabled = HOST_REGISTERED
        state%finalized = .FALSE.
    END SUBROUTINE add_hostcluster

    SUBROUTINE configure_hostcluster_kinematics(ntimes, t, xhost, yhost, zhost, vxhost, vyhost, vzhost)
        INTEGER, INTENT(IN) :: ntimes
        REAL*8, INTENT(IN), DIMENSION(ntimes) :: t, xhost, yhost, zhost, vxhost, vyhost, vzhost

        IF (.NOT. HOST_REGISTERED) CALL hostcluster_add()
        CALL hostcluster_configure_kinematics(ntimes, t, xhost, yhost, zhost, vxhost, vyhost, vzhost)
        CALL clear_force_registry()
        state%host_enabled = HOST_REGISTERED
        state%finalized = .FALSE.
        N_HOST_ORBIT_TIME_STAMPS = ntimes
    END SUBROUTINE configure_hostcluster_kinematics

    SUBROUTINE configure_hostcluster_structure(model_name, params, nparams)
        CHARACTER(LEN=*), INTENT(IN) :: model_name
        INTEGER, INTENT(IN) :: nparams
        REAL*8, INTENT(IN), DIMENSION(nparams) :: params

        IF (.NOT. HOST_REGISTERED) CALL hostcluster_add()
        CALL hostcluster_configure_model(model_name, params, nparams)
        CALL clear_force_registry()
        state%host_enabled = HOST_REGISTERED
        state%finalized = .FALSE.
        N_HOST_STRUCTURE_PARAMETERS = nparams
    END SUBROUTINE configure_hostcluster_structure

    SUBROUTINE finalize_hostcluster()
        CALL hostcluster_finalize()
        CALL clear_force_registry()
        state%host_enabled = HOST_REGISTERED
        state%finalized = .FALSE.
    END SUBROUTINE finalize_hostcluster    

    SUBROUTINE get_hostcluster_kinematics(ntimes, t, xhost, yhost, zhost, vxhost, vyhost, vzhost)
        INTEGER, INTENT(IN) :: ntimes
        REAL*8, INTENT(OUT), DIMENSION(ntimes) :: t, xhost, yhost, zhost, vxhost, vyhost, vzhost

        CALL hostcluster_get_kinematics(ntimes, t, xhost, yhost, zhost, vxhost, vyhost, vzhost)
    END SUBROUTINE get_hostcluster_kinematics

    SUBROUTINE get_hostcluster_structure(n_params,model_name,constant_params)
        INTEGER, INTENT(IN)                         :: n_params
        CHARACTER(LEN=64), INTENT(OUT)              :: model_name 
        REAL*8, INTENT(OUT), DIMENSION(n_params)    :: constant_params   
        CALL hostcluster_get_structure(n_params, model_name, constant_params)
    END SUBROUTINE get_hostcluster_structure 

    SUBROUTINE configure_hostcluster_structure_parameter_table(index,ntimes,timestamps_parameter,values)
        INTEGER, INTENT(IN) :: index, ntimes
        REAL*8, INTENT(IN), DIMENSION(ntimes) :: timestamps_parameter,values
        CALL hostcluster_configure_structure_parameter_table(index,ntimes,timestamps_parameter,values)
    END SUBROUTINE configure_hostcluster_structure_parameter_table

    subroutine get_hostcluster_ionization_state(NP, ionized_particles, ionization_time)
        INTEGER, INTENT(IN) :: NP
        LOGICAL, DIMENSION(NP), INTENT(OUT) :: ionized_particles
        REAL*8, DIMENSION(NP), INTENT(OUT) :: ionization_time
        IF (Nparticles /= NP) then 
            print*, "WARNING: get_hostcluster_ionization_state, nparticles mismatch"
        END IF 
        CALL hostcluster_get_ionization_state(NP, ionized_particles, ionization_time)
        
    end subroutine

    !!!! OUTPUTS
    SUBROUTINE initwritesnapshots(nskip, directory, basename)
        INTEGER, INTENT(IN) :: nskip
        CHARACTER(LEN=*), INTENT(IN) :: directory, basename
        LOGICAL :: dir_exists

        INQUIRE(FILE=directory, EXIST = dir_exists)

        if (.NOT.dir_exists) then             
            print*, "WANRING: initwritesnapshots. The directory", directory, "doesn't exist "
            print*, "   Create the directory manually before executing the code"
            print*, "           --- ABORTING ---"
            state%DONTRUN = .TRUE.
            RETURN
        end if 

        directory_snapshots = TRIM(directory)
        basename_snapshots = TRIM(basename)
        nskip_writesnapshots = nskip
        state%writesnapshots = .TRUE.
    END SUBROUTINE initwritesnapshots

    SUBROUTINE initwriteorbits(nskip, directory, basename)
        INTEGER, INTENT(IN) :: nskip
        CHARACTER(LEN=*), INTENT(IN) :: directory, basename

        directory_orbits = TRIM(directory)
        basename_orbits = TRIM(basename)
        nskip_writeorbits = nskip
        state%writeorbits = .TRUE.
    END SUBROUTINE initwriteorbits
        
    SUBROUTINE write_snapshot_file(istep)
        INTEGER, INTENT(IN) :: istep
        CHARACTER(LEN=600) :: filepath
        CHARACTER(LEN=20) :: step_str
        INTEGER :: iunit
        REAL*4, DIMENSION(:), ALLOCATABLE :: x_sp, y_sp, z_sp, vx_sp, vy_sp, vz_sp
        REAL*4 :: time_sp

        ! Construct filepath: directory/basename_ISTEP.bin
        WRITE(step_str, '(I0.6)') istep
        filepath = TRIM(directory_snapshots) // '/' // TRIM(basename_snapshots) // '_' // TRIM(ADJUSTL(step_str)) // '.bin'

        ! Convert to single precision for compact storage
        ALLOCATE(x_sp(Nparticles), y_sp(Nparticles), z_sp(Nparticles))
        ALLOCATE(vx_sp(Nparticles), vy_sp(Nparticles), vz_sp(Nparticles))
        x_sp = REAL(x(1:Nparticles), kind=4)
        y_sp = REAL(y(1:Nparticles), kind=4)
        z_sp = REAL(z(1:Nparticles), kind=4)
        vx_sp = REAL(vx(1:Nparticles), kind=4)
        vy_sp = REAL(vy(1:Nparticles), kind=4)
        vz_sp = REAL(vz(1:Nparticles), kind=4)
        time_sp = REAL(timestamps(current_step), kind=4)

        ! Open file for writing (unformatted binary, single precision)
        OPEN(NEWUNIT=iunit, FILE=TRIM(filepath), STATUS='REPLACE', ACTION='WRITE', FORM='UNFORMATTED')

        ! Write header: nparticles, time
        WRITE(iunit) Nparticles, time_sp

        ! Write phase-space data: x, y, z, vx, vy, vz (single precision)
        WRITE(iunit) x_sp
        WRITE(iunit) y_sp
        WRITE(iunit) z_sp
        WRITE(iunit) vx_sp
        WRITE(iunit) vy_sp
        WRITE(iunit) vz_sp

        CLOSE(iunit)
        DEALLOCATE(x_sp, y_sp, z_sp, vx_sp, vy_sp, vz_sp)

    END SUBROUTINE write_snapshot_file

    SUBROUTINE open_orbit_files()
        ! Open all orbit files once before integration loop (one per particle)
        ! No memory limit: all Nparticles get their own file
        ! Files remain open for the entire integration for fast writing
        CHARACTER(LEN=600) :: filepath
        CHARACTER(LEN=20) :: pid_str
        INTEGER :: ipart
        
        IF (ALLOCATED(orbit_file_units)) DEALLOCATE(orbit_file_units)
        ALLOCATE(orbit_file_units(Nparticles))
        orbit_file_units = 0
        
        DO ipart = 1, Nparticles
            WRITE(pid_str, '(I0.6)') ipart
            filepath = TRIM(directory_orbits) // '/' // TRIM(basename_orbits) // '_particle_' // TRIM(ADJUSTL(pid_str)) // '.bin'
            
            ! Open file for unformatted sequential write
            OPEN(NEWUNIT=orbit_file_units(ipart), FILE=TRIM(filepath), &
                 STATUS='UNKNOWN', ACTION='WRITE', FORM='UNFORMATTED')
        END DO
    
    END SUBROUTINE open_orbit_files
    
    SUBROUTINE close_orbit_files()
        ! Close all orbit files at the end of integration
        INTEGER :: ipart
        
        IF (ALLOCATED(orbit_file_units)) THEN
            DO ipart = 1, SIZE(orbit_file_units)
                IF (orbit_file_units(ipart) > 0) THEN
                    CLOSE(orbit_file_units(ipart))
                    orbit_file_units(ipart) = 0
                END IF
            END DO
        END IF
    
    END SUBROUTINE close_orbit_files

    SUBROUTINE write_orbit_records()
        ! Write current state to already-open per-particle orbit files
        ! One file per particle (unlimited, no memory constraint); each call appends one record
        INTEGER :: ipart
        REAL*4 :: time_sp, x_sp, y_sp, z_sp, vx_sp, vy_sp, vz_sp

        IF (.NOT. ALLOCATED(orbit_file_units)) RETURN
        
        DO ipart = 1, Nparticles
            ! Convert current state to single precision
            time_sp = REAL(timestamps(current_step), kind=4)
            x_sp = REAL(x(ipart), kind=4)
            y_sp = REAL(y(ipart), kind=4)
            z_sp = REAL(z(ipart), kind=4)
            vx_sp = REAL(vx(ipart), kind=4)
            vy_sp = REAL(vy(ipart), kind=4)
            vz_sp = REAL(vz(ipart), kind=4)

            ! Write single record to already-open unit: time, x, y, z, vx, vy, vz
            WRITE(orbit_file_units(ipart)) time_sp, x_sp, y_sp, z_sp, vx_sp, vy_sp, vz_sp
            FLUSH(orbit_file_units(ipart))
        END DO

    END SUBROUTINE write_orbit_records

    ! COMPUTATION AND PREPARATIONS
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
        IF (allocated(orbits_timestamps)) DEALLOCATE(orbits_timestamps)
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
        allocate(orbits_timestamps(NSAVED_ORBITS))
        state%orbits_allocated=.TRUE.

    END SUBROUTINE allocate_orbits

    ! THE NUMERICAL SCHEMES
    SUBROUTINE leapfrog()
        REAL*8, DIMENSION(Nparticles) :: xtemp, ytemp, ztemp
        REAL*8, DIMENSION(Nparticles) :: fx, fy, fz
        REAL*8 :: dt_step

        IF (current_step > nsteps) RETURN

        currenttime = timestamps(current_step)
        dt_step = timestamps(current_step + 1) - timestamps(current_step)

        xtemp = x + 0.5D0 * dt_step * vx
        ytemp = y + 0.5D0 * dt_step * vy
        ztemp = z + 0.5D0 * dt_step * vz

        CALL evaluate_total_force(currenttime + 0.5D0*dt_step, Nparticles, xtemp, ytemp, ztemp, fx, fy, fz)

        vx = vx + dt_step * fx
        vy = vy + dt_step * fy
        vz = vz + dt_step * fz

        x = xtemp + 0.5D0 * dt_step * vx
        y = ytemp + 0.5D0 * dt_step * vy
        z = ztemp + 0.5D0 * dt_step * vz

        current_step = current_step + 1
        currenttime = timestamps(current_step)
    END SUBROUTINE leapfrog

    SUBROUTINE forest_ruth()
        REAL*8, DIMENSION(Nparticles) :: fx, fy, fz
        REAL*8 :: dt_step

        IF (current_step > nsteps) RETURN

        dt_step = timestamps(current_step + 1) - timestamps(current_step)
        
        currenttime = timestamps(current_step)

        ! DRIFT
        x = x + yoshida_c1 * dt_step * vx
        y = y + yoshida_c1 * dt_step * vy
        z = z + yoshida_c1 * dt_step * vz
        ! KICK 
        currenttime = currenttime + yoshida_c1*dt_step 
        CALL evaluate_total_force(currenttime, Nparticles, x, y, z, fx, fy, fz)
        vx = vx + yoshida_d1*fx*dt_step
        vy = vy + yoshida_d1*fy*dt_step
        vz = vz + yoshida_d1*fz*dt_step
        ! DRIFT
        x = x + yoshida_c2 * dt_step * vx
        y = y + yoshida_c2 * dt_step * vy
        z = z + yoshida_c2 * dt_step * vz
        ! KICK 
        currenttime = currenttime + yoshida_c2*dt_step 
        CALL evaluate_total_force(currenttime, Nparticles, x, y, z, fx, fy, fz)
        vx = vx + yoshida_d2*fx*dt_step
        vy = vy + yoshida_d2*fy*dt_step
        vz = vz + yoshida_d2*fz*dt_step
        ! DRIFT
        x = x + yoshida_c3 * dt_step * vx
        y = y + yoshida_c3 * dt_step * vy
        z = z + yoshida_c3 * dt_step * vz
        ! KICK 
        currenttime = currenttime + yoshida_c3*dt_step 
        CALL evaluate_total_force(currenttime, Nparticles, x, y, z, fx, fy, fz)
        vx = vx + yoshida_d3*fx*dt_step
        vy = vy + yoshida_d3*fy*dt_step
        vz = vz + yoshida_d3*fz*dt_step
        ! DRIFT
        x = x + yoshida_c4 * dt_step * vx
        y = y + yoshida_c4 * dt_step * vy
        z = z + yoshida_c4 * dt_step * vz
        ! KICK 
        currenttime = currenttime + yoshida_c4*dt_step 
        CALL evaluate_total_force(currenttime, Nparticles, x, y, z, fx, fy, fz)
        vx = vx + yoshida_d4*fx*dt_step
        vy = vy + yoshida_d4*fy*dt_step
        vz = vz + yoshida_d4*fz*dt_step

        current_step = current_step + 1 
        currenttime = timestamps(current_step)


    END SUBROUTINE forest_ruth

    SUBROUTINE compute_yoshida_coefficients()
        yoshida_w = sqrt(2.0D0**(1.0D0/3.0D0) + 2.0D0**(-1.0D0/3.0D0) -1.0D0 )/6.0D0 ! D0 is for double precision
        yoshida_c1 =  yoshida_w + 0.5D0
        yoshida_c2 = -yoshida_w
        yoshida_c3 = -yoshida_w
        yoshida_c4 =  yoshida_w + 0.5D0    
        yoshida_d1 =  2.0D0*yoshida_w+1.0D0
        yoshida_d2 = -4.0D0*yoshida_w-1.0D0
        yoshida_d3 =  2.0D0*yoshida_w+1.0D0
        yoshida_d4 =  0.0D0        
    END SUBROUTINE compute_yoshida_coefficients

    ! FORCE ORCHESTRATION aross multiple modules
    SUBROUTINE clear_force_registry()
        INTEGER :: i

        nactive_force_providers = 0
        DO i = 1, MAX_ACTIVE_FORCE_PROVIDERS
            active_force_providers(i)%name = ""
            NULLIFY(active_force_providers(i)%eval)
        END DO
    END SUBROUTINE clear_force_registry

    SUBROUTINE register_force_provider(name, eval_proc)
        CHARACTER(LEN=*), INTENT(IN) :: name
        PROCEDURE(force_provider_interface) :: eval_proc

        IF (nactive_force_providers >= MAX_ACTIVE_FORCE_PROVIDERS) THEN
            PRINT*, "WARNING: register_force_provider: maximum providers reached"
            RETURN
        END IF

        nactive_force_providers = nactive_force_providers + 1
        active_force_providers(nactive_force_providers)%name = TRIM(name)
        active_force_providers(nactive_force_providers)%eval => eval_proc
    END SUBROUTINE register_force_provider

    SUBROUTINE rebuild_force_registry()
        CALL clear_force_registry()

        IF (GRAVITY_FINALIZED .AND. GRAVITY_NCOMP > 0) THEN
            CALL register_force_provider("gravity", gravity_force_provider)
        END IF

        IF (HOST_REGISTERED .AND. HOST_FINALIZED) THEN
            CALL register_force_provider("hostcluster", hostcluster_force_provider)
        END IF

        IF (FLYINGSPHERES_REGISTERED .AND. FLYINGSPHERES_FINALIZED) THEN 
            CALL register_force_provider("flyingspheres", flyingspheres_force_provider)
        END IF 
    END SUBROUTINE rebuild_force_registry

    SUBROUTINE gravity_force_provider(t, n, xin, yin, zin, ax, ay, az)
        REAL*8, INTENT(IN) :: t
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: xin, yin, zin
        REAL*8, INTENT(OUT), DIMENSION(n) :: ax, ay, az

        IF (t /= t) THEN
            ax = 0.0D0
            ay = 0.0D0
            az = 0.0D0
            RETURN
        END IF

        CALL gravity_force(n, xin, yin, zin, ax, ay, az)
    END SUBROUTINE gravity_force_provider

    SUBROUTINE hostcluster_force_provider(t, n, xin, yin, zin, ax, ay, az)
        ! this ensures that the host state is updated each time the force is called
        REAL*8, INTENT(IN) :: t
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: xin, yin, zin
        REAL*8, INTENT(OUT), DIMENSION(n) :: ax, ay, az

        CALL hostcluster_update_state(t)
        CALL hostcluster_eval_force(n, xin, yin, zin, ax, ay, az)
    END SUBROUTINE hostcluster_force_provider

    SUBROUTINE flyingspheres_force_provider(t, n, xin, yin, zin, ax, ay, az)
        REAL*8, INTENT(IN) :: t
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: xin, yin, zin
        REAL*8, INTENT(OUT), DIMENSION(n) :: ax, ay, az
        call system_clock(c_flyingspheres_start)
        CALL flyingspheres_update_state(t)   
        call flyingspheres_force(n,xin,yin,zin,ax,ay,az) 
        call system_clock(c_flyingspheres_end)
        c_flyingspheres_accum = c_flyingspheres_accum + (c_flyingspheres_end-c_flyingspheres_start)

    END SUBROUTINE flyingspheres_force_provider

    SUBROUTINE evaluate_total_force(t, n, xin, yin, zin, ax, ay, az)
        REAL*8, INTENT(IN) :: t
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: xin, yin, zin
        REAL*8, INTENT(OUT), DIMENSION(n) :: ax, ay, az
        REAL*8, DIMENSION(n) :: ax_tmp, ay_tmp, az_tmp
        INTEGER :: i

        ax = 0.0D0
        ay = 0.0D0
        az = 0.0D0

        DO i = 1, nactive_force_providers
            IF (.NOT. ASSOCIATED(active_force_providers(i)%eval)) CYCLE
            CALL active_force_providers(i)%eval(t, n, xin, yin, zin, ax_tmp, ay_tmp, az_tmp)
            ax = ax + ax_tmp
            ay = ay + ay_tmp
            az = az + az_tmp
        END DO
    END SUBROUTINE evaluate_total_force

END MODULE simulator

