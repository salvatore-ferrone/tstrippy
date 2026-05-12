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
        LOGICAL :: writesnapshots = .FALSE.
        LOGICAL :: writeorbits = .FALSE.
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

    REAL*8, PARAMETER :: DEFAULT_ORBIT_RAM_LIMIT_MB = 1024.0D0
    REAL*8, PARAMETER :: DEFAULT_MAX_RAM_MB = 1024.0D0

    ! for saving some trajectories
    REAL*8, PUBLIC :: orbit_ram_limit_MB = DEFAULT_ORBIT_RAM_LIMIT_MB
    INTEGER, PRIVATE :: n_particles_orbit = 1 
    INTEGER, PRIVATE :: nskip_orbit_timestamps = 1 
    INTEGER, PRIVATE :: nvars_orbits = 7 ! time and phase space 
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: orbits
    
    ! persistent file handles for direct orbit writing (no per-step open/close)
    INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: orbit_file_units

    ! i/o
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

    CONTAINS 
    !!!!! CALLS WHERE THE USER INTERFACES WITH THE MODULE
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
        IF (ALLOCATED(orbits)) DEALLOCATE(orbits)
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

        IF (ALLOCATED(orbits)) DEALLOCATE(orbits)
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
        state%orbits_allocated = .FALSE.
        state%timestamps_from_user = .TRUE.
        state%finalized = .FALSE.
    END SUBROUTINE settimestamps

    SUBROUTINE setbackwardorbit()
        state%backward_orbit = .TRUE.
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
        state%orbits_allocated = .FALSE.
        state%finalized = .FALSE.
    END SUBROUTINE trim_orbits

    !!!! OUTPUTS
    SUBROUTINE initwritesnapshots(nskip, directory, basename)
        INTEGER, INTENT(IN) :: nskip
        CHARACTER(LEN=*), INTENT(IN) :: directory, basename

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
        IF (ALLOCATED(orbits)) DEALLOCATE(orbits)

        ! Scheme pointer and name
        NULLIFY(scheme)
        scheme_name = ""

        ! Counters
        Nparticles   = 0
        nsteps       = 0
        current_step = 1
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

        ! Reset all state flags to defaults
        state = state_t()

    END SUBROUTINE CLEAR

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

        nsteps = SIZE(timestamps) - 1
        current_step = 1

        if (should_return) return

        if (.not.state%orbits_allocated) THEN
            CALL allocate_orbits(nskip_orbit_timestamps)
            IF (.NOT. state%orbits_allocated) THEN
                PRINT*, "ERROR: finalize failed because orbit allocation did not complete"
                RETURN
            END IF
        END IF

        state%scheme_set = .TRUE.
        state%finalized = .TRUE.

    END subroutine finalize

    SUBROUTINE run()
        INTEGER :: istep, iorbit

        CALL finalize()

        IF (.NOT. state%finalized) THEN
            PRINT*, "ERROR: finalize failed. Cannot Run"
            RETURN
        END IF
        
        ! Open all orbit files before loop if writing orbits (one per particle, unlimited)
        IF (state%writeorbits .AND. Nparticles > 0) THEN
            CALL open_orbit_files()
            ! Write initial conditions immediately
            CALL write_orbit_records()
        END IF

        ! Save initial conditions as the first orbit snapshot
        iorbit = 1
        IF (n_particles_orbit > 0) THEN
            orbits(iorbit, 1, :) = timestamps(1)
            orbits(iorbit, 2, :) = x(1:n_particles_orbit)
            orbits(iorbit, 3, :) = y(1:n_particles_orbit)
            orbits(iorbit, 4, :) = z(1:n_particles_orbit)
            orbits(iorbit, 5, :) = vx(1:n_particles_orbit)
            orbits(iorbit, 6, :) = vy(1:n_particles_orbit)
            orbits(iorbit, 7, :) = vz(1:n_particles_orbit)
            iorbit = iorbit + 1
        END IF

        ! Main integration loop
        DO istep = 1, nsteps
            CALL scheme()  ! advances positions/velocities and increments current_step

            ! Save orbit snapshot every nskip_orbit_timestamps steps
            IF (n_particles_orbit > 0) THEN
                IF (MOD(istep, nskip_orbit_timestamps) == 0) THEN
                    orbits(iorbit, 1, :) = timestamps(current_step)
                    orbits(iorbit, 2, :) = x(1:n_particles_orbit)
                    orbits(iorbit, 3, :) = y(1:n_particles_orbit)
                    orbits(iorbit, 4, :) = z(1:n_particles_orbit)
                    orbits(iorbit, 5, :) = vx(1:n_particles_orbit)
                    orbits(iorbit, 6, :) = vy(1:n_particles_orbit)
                    orbits(iorbit, 7, :) = vz(1:n_particles_orbit)
                    iorbit = iorbit + 1
                END IF
            END IF

            ! Write snapshot if enabled
            IF (state%writesnapshots) THEN
                IF (MOD(istep, nskip_writesnapshots) == 0) THEN
                    CALL write_snapshot_file(istep)
                END IF
            END IF

            ! Write orbits directly (per-particle files, no memory limit)
            IF (state%writeorbits .AND. ALLOCATED(orbit_file_units)) THEN
                IF (MOD(istep, nskip_writeorbits) == 0) THEN
                    CALL write_orbit_records()
                END IF
            END IF
        END DO
        
        ! Close all orbit files at end
        IF (ALLOCATED(orbit_file_units)) THEN
            CALL close_orbit_files()
        END IF

    END SUBROUTINE run

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
