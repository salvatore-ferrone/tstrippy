MODULE simulator 
    ! this integrator needs to contain the current positions
    ! it needs to be able to apply any force that I want at any time
    ! it needs to be able to integrate the positions and velocities
    use gravity, ONLY: pot_clear_basis                                   => clearsphericalharmonicbasis,                 &
                          pot_init_basis                                    => initsphericalharmonicbasis,                  &
                          pot_clear_composite_gravity                       => clearcompositegravity,                       &
                          pot_init_composite_gravity                        => initcompositegravity,                        &
                          pot_add_composite_sphericalharmonic_exponential_oblate => addcompositesphericalharmonicexponentialoblate, &
                          pot_add_composite_sphericalharmonic_ibata2024_halo     => addcompositesphericalharmonicibata2024halo,     &
                          pot_add_composite_disk_bessel_exponential_disk         => addcompositediskbesselexponentialdisk,           &
                          pot_finalize_composite_gravity                    => finalizecompositegravity,                    &
                          grav_cleargravity                                 => cleargravity,                                &
                          grav_setgravityconstant                           => setgravityconstant,                          &
                          grav_addgravitycomponent                          => addgravitycomponent,                         &
                          grav_finalizegravity                              => finalizegravity,                             &
                          grav_evaluategravityforces                        => evaluategravityforces,                       &
                          grav_evaluategravitypotential                     => evaluategravitypotential,                    &
                          grav_gravity_g                                    => GRAVITY_G,                                   &
                          grav_gravity_finalized                            => GRAVITY_FINALIZED,                           &
                          pot_nbody_plummers                                => NBODYPLUMMERS
    use perturbers, ONLY: pert_init          => perturberinitialization, &
                          pert_deallocate    => perturberdeallocation,   &
                          pert_find_time_idx => findperturbertimeindex,  &
                          pert_compute_force => computeforcebyperturbers
    use hostperturber, ONLY: host_init_kinematics_mod => host_init_kinematics, &
                             host_init_mass_mod       => host_init_mass,       &
                             host_init_radius_mod     => host_init_radius,     &
                             host_deallocate_mod      => hostdeallocation,     &
                             host_find_time_idx       => findhosttimeindex,    &
                             host_compute_force       => computeforcebyhosts,  &
                             host_vx_current          => vxhostcurrent,        &
                             host_vy_current          => vyhostcurrent,        &
                             host_vz_current          => vzhostcurrent
    use galacticbar, ONLY: bar_init          => galacticbarinitialization, &
                          bar_deallocate    => bardeallocation,          &
                          bar_update        => updatebarorientation,     &
                          bar_force_eval    => barforce,                 &
                          bar_angle         => theta
    IMPLICIT NONE
    PRIVATE 
    ! DECLARE SUBROUTINES
    PUBLIC :: setgravityconstant, cleargravitycomponents, addgravitycomponent, finalizegravity
    PUBLIC :: finalizesimulator, evaluatepotentials
    PUBLIC :: setstaticgalaxy, setintegrationparameters, setinitialkinematics
    PUBLIC :: setdebugaccelerations, setdebugbarorientation, setbackwardorbit
    PUBLIC :: inithostkinematics, inithostmass, inithostradius
    PUBLIC :: initnbodysystem,initgalacticbar,initperturbers
    PUBLIC :: leapfrogintime, leapfrogtofinalpositions
    PUBLIC :: ruthforestintime
    PUBLIC :: setintegratormethod, runsimulation
    PUBLIC :: HIT
    PUBLIC :: assert_gravitational_constant_initialized
    PUBLIC :: initwriteparticleorbits, writeparticleorbits
    PUBLIC :: initwritestream, writestream
    PUBLIC :: initwritesnapshot, writesnapshot
    PUBLIC :: deallocate
    PUBLIC :: clearsphericalharmonicbasis, initsphericalharmonicbasis
    PUBLIC :: clearcompositegravity, initcompositegravity
    PUBLIC :: addcompositesphericalharmonicexponentialoblate, addcompositesphericalharmonicibata2024halo
    PUBLIC :: addcompositediskbesselexponentialdisk
    PUBLIC :: finalizecompositegravity
    ! DECIDE WHICH PHYSICS TO INCLUDE
    LOGICAL, PUBLIC :: DONBODY = .FALSE.
    LOGICAL, PUBLIC :: DOPERTURBERS = .FALSE.
    LOGICAL, PUBLIC :: DOHOSTPERTURBER = .FALSE.
    LOGICAL, PUBLIC :: DOGALACTICBAR = .FALSE.
    LOGICAL, PUBLIC :: DOBACKWARDORBIT = .FALSE.
    ! Variables to keep track of the physics that has been set
    LOGICAL, PUBLIC :: INITIALKINEMATICSSET = .FALSE.
    LOGICAL, PUBLIC :: INTEGRATIONPARAMETERSSET = .FALSE.
    ! DECIDE IF WE WILL BE SAVING WHOLE ORBITS OR SNAPSHOTS
    LOGICAL, PUBLIC :: DOWRITEORBITS = .FALSE.
    LOGICAL, PUBLIC :: DOWRITESTREAM = .FALSE.
    ! DEBUGGING VARIABLES 
    LOGICAL, PUBLIC :: DEBUGACCELERATIONS = .FALSE. ! save the accelerations for debugging
    REAL*8, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: aSG,aHP,aP,aNBODY,aBAR,aTOTAL
    LOGICAL, PUBLIC :: DEBUGBARORIENTATION = .FALSE.
    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: bartheta
    ! DECLARE MODULE WIDE VARIABLES
    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: timestamps
    REAL*8,DIMENSION(:),ALLOCATABLE,PUBLIC :: xf,yf,zf,vxf,vyf,vzf,tesc,nbodyparams
    REAL*8, PUBLIC :: currenttime,dt
    INTEGER, PUBLIC :: ntimesteps,ntimepoints,nparticles,nwriteskip
    INTEGER, PUBLIC :: INTEGRATIONMETHOD = 0 ! 0=leapfrog, 1=forest_ruth
    INTEGER, PUBLIC :: FILEUNITBASE 
    CHARACTER*500, PUBLIC :: outname,outdir,streamdir,streamname
    !! THE ACCELEARTIONS ARE PUBLIC SO THAT THEY CAN BE ACCESSED BY THE DEBUGGING SUBROUTINES
    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: axSG,aySG,azSG
    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: axHP,ayHP,azHP
    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: axP,ayP,azP
    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: axNBODY,ayNBODY,azNBODY
    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: axBAR,ayBAR,azBAR
    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: phiSG,phiHP,phiP,phiBAR
    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: phiNBODY
    REAL*8, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: phiTensor
    contains 
    SUBROUTINE setgravityconstant(Gin)
        REAL*8, INTENT(IN) :: Gin
        CALL grav_setgravityconstant(Gin)
    END SUBROUTINE setgravityconstant

    SUBROUTINE cleargravitycomponents()
        CALL grav_cleargravity()
    END SUBROUTINE cleargravitycomponents

    SUBROUTINE addgravitycomponent(modelname, params)
        CHARACTER(LEN=*), INTENT(IN) :: modelname
        REAL*8, DIMENSION(:), INTENT(IN) :: params
        CALL grav_addgravitycomponent(modelname, params, SIZE(params))
    END SUBROUTINE addgravitycomponent

    SUBROUTINE finalizegravity()
        CALL grav_finalizegravity()
    END SUBROUTINE finalizegravity

    SUBROUTINE setstaticgalaxy(milkywaypotentialname,mwparams)
        ! Legacy wrapper. Expects mwparams = [G, model_params...].
        CHARACTER*100, INTENT(IN) :: milkywaypotentialname
        REAL*8, DIMENSION(:), INTENT(IN) :: mwparams

        IF (SIZE(mwparams) < 2) THEN
            PRINT*, "ERROR: setstaticgalaxy requires [G, model_params...]"
            STOP
        END IF

        CALL cleargravitycomponents()
        CALL setgravityconstant(mwparams(1))
        CALL addgravitycomponent(TRIM(milkywaypotentialname), mwparams(2:))
        CALL finalizegravity()
    END SUBROUTINE setstaticgalaxy

    ! SUBROUTINES FOR MANUALLY SETTING THE BASIS EXPANSION PARAMS
    SUBROUTINE clearsphericalharmonicbasis()
        CALL pot_clear_basis()
    END SUBROUTINE clearsphericalharmonicbasis

    SUBROUTINE initsphericalharmonicbasis(lmax, nr, r_grid)
        INTEGER, INTENT(IN) :: lmax, nr
        REAL*8, DIMENSION(nr), INTENT(IN) :: r_grid
        CALL pot_init_basis(lmax, nr, r_grid)
    END SUBROUTINE initsphericalharmonicbasis

    SUBROUTINE clearcompositegravity()
        CALL pot_clear_composite_gravity()
    END SUBROUTINE clearcompositegravity

    SUBROUTINE initcompositegravity(lmax, nr, r_grid, ncomp)
        INTEGER, INTENT(IN) :: lmax, nr, ncomp
        REAL*8, DIMENSION(nr), INTENT(IN) :: r_grid
        CALL pot_init_composite_gravity(lmax, nr, r_grid, ncomp)
    END SUBROUTINE initcompositegravity

    SUBROUTINE addcompositesphericalharmonicexponentialoblate(component_index, rho0, s0, q)
        INTEGER, INTENT(IN) :: component_index
        REAL*8, INTENT(IN) :: rho0, s0, q
        CALL pot_add_composite_sphericalharmonic_exponential_oblate(component_index, rho0, s0, q)
    END SUBROUTINE addcompositesphericalharmonicexponentialoblate

    SUBROUTINE addcompositesphericalharmonicibata2024halo(component_index, rho0, r0, rt, q, gamma, beta)
        INTEGER, INTENT(IN) :: component_index
        REAL*8, INTENT(IN) :: rho0, r0, rt, q, gamma, beta
        CALL pot_add_composite_sphericalharmonic_ibata2024_halo(component_index, rho0, r0, rt, q, gamma, beta)
    END SUBROUTINE addcompositesphericalharmonicibata2024halo

    SUBROUTINE addcompositediskbesselexponentialdisk(component_index, sigma0, hR, hZ)
        INTEGER, INTENT(IN) :: component_index
        REAL*8, INTENT(IN) :: sigma0, hR, hZ
        CALL pot_add_composite_disk_bessel_exponential_disk(component_index, sigma0, hR, hZ)
    END SUBROUTINE addcompositediskbesselexponentialdisk

    SUBROUTINE finalizecompositegravity()
        CALL pot_finalize_composite_gravity()
    END SUBROUTINE finalizecompositegravity

    SUBROUTINE assert_gravitational_constant_initialized()
        if (.NOT. grav_gravity_finalized) then
            print*, "E300: gravity not finalized before integration"
            print*, "      call finalizegravity() or finalizesimulator() first"
            stop
        END IF
    END SUBROUTINE assert_gravitational_constant_initialized

    SUBROUTINE setinitialkinematics(N,x,y,z,vx,vy,vz)
        ! set the initial kinematics of the particles
        INTEGER, intent(in) :: N
        REAL*8, DIMENSION(N), intent(in) :: x,y,z,vx,vy,vz
        allocate(xf(N),yf(N),zf(N),vxf(N),vyf(N),vzf(N),tesc(N))
        allocate(axSG(N),aySG(N),azSG(N))
        allocate(axHP(N),ayHP(N),azHP(N))
        allocate(axP(N),ayP(N),azP(N))
        allocate(axNBODY(N),ayNBODY(N),azNBODY(N))
        allocate(axBAR(N),ayBAR(N),azBAR(N))
        allocate(phiSG(N),phiHP(N),phiP(N),phiBAR(N),phiNBODY(N))
        allocate(phiTensor(N,N))
        nparticles=N
        xf = x
        yf = y
        zf = z
        vxf = vx
        vyf = vy
        vzf = vz
        tesc = -9990.0
        INITIALKINEMATICSSET = .TRUE.
    END SUBROUTINE setinitialkinematics

    SUBROUTINE setintegrationparameters(t0,dt0,nsteps)
        ! define the total integration time, the timestep, and the number of timesteps
        REAL*8, intent(in) :: t0,dt0
        INTEGER, intent(in) :: nsteps
        integer :: i
        currenttime = t0
        dt = dt0
        ntimesteps = nsteps
        ntimepoints = nsteps + 1
        allocate(timestamps(ntimepoints))
        timestamps(1) = t0
        DO i=2,ntimepoints
            timestamps(i) = timestamps(i-1) + dt
        END DO
        INTEGRATIONPARAMETERSSET = .TRUE.

        
    END SUBROUTINE setintegrationparameters

    SUBROUTINE setintegratormethod(methodname)
        ! Select integration backend for runsimulation().
        CHARACTER*100, INTENT(IN) :: methodname

        IF (TRIM(methodname) .EQ. "leapfrog") THEN
            INTEGRATIONMETHOD = 0
        ELSE IF (TRIM(methodname) .EQ. "forest_ruth") THEN
            INTEGRATIONMETHOD = 1
        ELSE
            PRINT*, "E301: invalid integration method"
            PRINT*, "      valid methods: leapfrog, forest_ruth"
            STOP
        END IF
    END SUBROUTINE setintegratormethod

    SUBROUTINE runsimulation()
        ! Unified run entrypoint (first refactor slice).
        ! Consumes configured module state and updates final phase-space arrays.
        IF (.NOT. grav_gravity_finalized) THEN
            PRINT*, "E302: gravitational field not finalized before runsimulation"
            PRINT*, "      call finalizegravity() or finalizesimulator() first"
            STOP
        END IF
        IF (.NOT. INITIALKINEMATICSSET) THEN
            PRINT*, "E303: setinitialkinematics must be called before runsimulation"
            STOP
        END IF
        IF (.NOT. INTEGRATIONPARAMETERSSET) THEN
            PRINT*, "E304: setintegrationparameters must be called before runsimulation"
            STOP
        END IF

        IF (INTEGRATIONMETHOD .EQ. 0) THEN
            CALL leapfrogtofinalpositions()
        ELSE IF (INTEGRATIONMETHOD .EQ. 1) THEN
            PRINT*, "E305: runsimulation forest_ruth backend not wired in this slice"
            PRINT*, "      use ruthforestintime directly for now"
            STOP
        ELSE
            PRINT*, "E306: unknown INTEGRATIONMETHOD state"
            STOP
        END IF
    END SUBROUTINE runsimulation
    
    SUBROUTINE setbackwardorbit()
        ! Changes the sign of the velocities and the timestamps
        ! the timestamps take the current time and subtract dt from it over NSTEPS
        integer :: i

        if (INITIALKINEMATICSSET.eqv..FALSE.) then
            print*, "ERROR: setinitialkinematics must be called before setbackwardorbit"
            stop
        end if
        vxf = -vxf
        vyf = -vyf
        vzf = -vzf
        ! reset the timestamps to go backward
        timestamps(1) = currenttime
        DO i=2,ntimepoints
            timestamps(i) = timestamps(i-1) - dt
        END DO
        DOBACKWARDORBIT = .TRUE.
    END SUBROUTINE setbackwardorbit

    SUBROUTINE setdebugaccelerations()
        if (INITIALKINEMATICSSET .eqv. .FALSE.) then
            print*, "ERROR: setinitialkinematics must be called before setdebugaccelerations"
            stop
        end if
        
        if (nparticles.ne.1) then
            print*, "ERROR: DEBUGACCELERATIONS only works for one particle"
            stop
        end if
        DEBUGACCELERATIONS = .TRUE.
        allocate(aSG(3,ntimepoints),aHP(3,ntimepoints),aP(3,ntimepoints))
        allocate(aNBODY(3,ntimepoints),aBAR(3,ntimepoints),aTOTAL(3,ntimepoints))
    END SUBROUTINE setdebugaccelerations

    SUBROUTINE setdebugbarorientation()
        if (INTEGRATIONPARAMETERSSET .eqv. .FALSE.) then
            print*, "ERROR: setintegrationparameters must be called before setdebugbarorientation"
            stop
        end if
        if (DOGALACTICBAR.eqv..FALSE.) then
            print*, "ERROR: initgalacticbar must be called before setdebugbarorientation"
            stop
        end if
        DEBUGBARORIENTATION = .TRUE.
        allocate(bartheta(ntimepoints))
    END SUBROUTINE setdebugbarorientation
    
    subroutine initnbodysystem(N,Gin,massesnbody,scaleradiinbody)
        ! initialize the nbody system
        ! meaning that everytime the system is evaluated, we also compute the Nbody forces
        INTEGER, intent(in) :: N
        REAL*8, DIMENSION(N), intent(in) :: massesnbody,scaleradiinbody
        REAL*8, intent(in) :: Gin
        DONBODY = .TRUE.
        allocate(nbodyparams(2*N+1))
        nbodyparams(1)=Gin
        nbodyparams(1+1:N+1)=massesnbody
        nbodyparams(N+1+1:2*N+1)=scaleradiinbody
    end subroutine initnbodysystem

    subroutine inithostmass(mass_model_name, params)
        ! sets the host mass evolution model. 
        character*100, intent(in) :: mass_model_name
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        INTEGER :: mass_model

        if (mass_model_name.eq."constant") then
            mass_model = 0
        else if (mass_model_name.eq."double_exponential") then
            mass_model = 1
        else
            print*, "ERROR: mass model not recognized in sethostmass"
            print*, "       available models are: constant, double_exponential"
            stop
        end if
        CALL host_init_mass_mod(mass_model, params)

    END SUBROUTINE inithostmass

    subroutine inithostradius(radius)
        ! sets the host radius
        REAL*8, intent(in) :: radius
        CALL host_init_radius_mod(radius)
    END SUBROUTINE inithostradius

    subroutine inithostkinematics(nhosttimepoints,timeH,xH,yH,zH,vxH,vyH,vzH)
        ! initialize the host perturber
        INTEGER, intent(in) :: nhosttimepoints
        real*8, intent(in), dimension(nhosttimepoints) :: timeH,xH,yH,zH,vxH,vyH,vzH
        DOHOSTPERTURBER = .TRUE.
        CALL host_init_kinematics_mod(nhosttimepoints,timeH,xH,yH,zH,vxH,vyH,vzH)

    end subroutine inithostkinematics

    SUBROUTINE initperturbers(tp,xp,yp,zp,masses,radii)
        real*8, intent(in) ,dimension(:) :: tp
        real*8, intent(in), dimension(:,:) :: xp,yp,zp
        real*8, intent(in), dimension(:) :: masses,radii
        DOPERTURBERS = .TRUE.
        CALL pert_init(SIZE(masses,1),SIZE(tp),tp,xp,yp,zp,masses,radii)
    END SUBROUTINE initperturbers

    SUBROUTINE initgalacticbar(barpotenname,barparams,barpoly)
        ! initialize the galactic bar
        character*100, INTENT(IN) :: barpotenname
        REAL*8, DIMENSION(:), intent(in) :: barparams,barpoly
        DOGALACTICBAR = .TRUE.
        CALL bar_init(barpotenname,barparams,barpoly)
    END SUBROUTINE initgalacticbar

    SUBROUTINE initwritestream(nskip,myoutname,myoutdir,memorybaseint)
        ! each time step gets its own file, contrary to writeparticleorbits
        ! the file will only have the positions and velocities
        ! the file will be written every NSKIP STEPS
        ! the files will be named outname-1.bin, outname-2.bin, etc.
        ! thie should be incompatible with writeparticleorbits
        ! I should block them from happening at the same time somehow
        INTEGER, intent(in) :: nskip
        CHARACTER*500, intent(in) :: myoutname
        CHARACTER*500, intent(in) :: myoutdir
        integer, intent(in),optional :: memorybaseint ! memory address to start the file units


        if (present(memorybaseint)) then
            FILEUNITBASE=memorybaseint 
        else
            FILEUNITBASE=10000
        end if

        DOWRITESTREAM=.TRUE.
        streamname=trim(myoutname)
        streamdir=trim(myoutdir)
        nwriteskip=nskip
    END SUBROUTINE initwritestream

    SUBROUTINE initwritesnapshot(nskip,myoutname,myoutdir,memorybaseint)
        ! New naming alias for stream-style snapshot writing.
        INTEGER, intent(in) :: nskip
        CHARACTER*500, intent(in) :: myoutname
        CHARACTER*500, intent(in) :: myoutdir
        integer, intent(in),optional :: memorybaseint

        CALL initwritestream(nskip,myoutname,myoutdir,memorybaseint)
    END SUBROUTINE initwritesnapshot

    SUBROUTINE writestream(index,N,x,y,z,vx,vy,vz)
        integer, intent(in) :: N
        REAL*8, DIMENSION(N), intent(in) :: x,y,z,vx,vy,vz
        integer::index
        character*500 :: filename
        ! make the file name
        WRITE(filename,'(A,A,A,I0,A)') &
            trim(streamdir),trim(streamname),'-',index,'.bin'
        ! open the file
        open(unit=FILEUNITBASE, file=filename, form="unformatted", status="replace")
        ! write the first line
        write(FILEUNITBASE) 6,N
        ! write the data
        write(FILEUNITBASE) SNGL(x),SNGL(y),SNGL(z),SNGL(vx),SNGL(vy),SNGL(vz)
        ! close the file
        close(FILEUNITBASE)
    END SUBROUTINE writestream

    SUBROUTINE writesnapshot(index,N,x,y,z,vx,vy,vz)
        ! New naming alias for stream-style snapshot writing.
        integer, intent(in) :: N
        REAL*8, DIMENSION(N), intent(in) :: x,y,z,vx,vy,vz
        integer::index

        CALL writestream(index,N,x,y,z,vx,vy,vz)
    END SUBROUTINE writesnapshot
    
    SUBROUTINE initwriteparticleorbits(nskip,myoutname,myoutdir,memorybaseint)
        ! INITIALIZE THE WRITING OF THE PARTICLE ORBITS
        ! EACH PARTICLE WILL HAVE ITS OWN FILE
        ! THE FILE WILL CONTAIN THE T,X,Y,Z,VX,VY,VZ
        ! THE FILE WILL BE WRITTEN EVERY NSKIP STEPS
        ! THE FILES WILL BE NAMED outname-1.bin, outname-2.bin, etc.
        ! THE FIRST LINE WILL BE THE (NSTEPS,8)
        INTEGER, intent(in) :: nskip
        integer, intent(in),optional :: memorybaseint ! memory address to start the file units
        CHARACTER*100, intent(in) :: myoutname
        CHARACTER*100, intent(in) :: myoutdir
        character*100 :: filename
        integer::i
        integer::nout,modulus
        if (present(memorybaseint)) then
            FILEUNITBASE=memorybaseint 
        else
            FILEUNITBASE=1000
        end if
        outname=trim(myoutname)
        outdir=trim(myoutdir)
        DOWRITEORBITS=.TRUE.
        nwriteskip=nskip
        modulus = mod(ntimepoints,nskip)
        if ((modulus).eq.(0)) then
            nout = ntimepoints/nskip
        ELSE
            ! if nskip is not a factor of the number of timesteps, 
            !   then we need to add one to include the last timestep
            nout = ntimepoints/nskip + 1 
        END if

        DO i=1,nparticles
            ! make the file name 
            WRITE(filename,'(A,A,A,I0,A)') trim(outdir),trim(outname),'-',i,'.dat'
            ! WRITE(filename,'(A,A,A,I0,A)') trim(outdir),trim(outname),'-',i,'.bin'
            ! open the file
            open(unit=FILEUNITBASE+i, file=filename, form="formatted", status="replace")
            ! write the first line
            ! write(FILEUNITBASE+i) nout,7
            ! close the file
            ! close(FILEUNITBASE+i)
        END DO

    END subroutine initwriteparticleorbits

    SUBROUTINE writeparticleorbits(myt,N,x,y,z,vx,vy,vz)
        integer, intent(in) :: N
        real*8, intent(in) :: myt
        REAL*8, DIMENSION(N), intent(in) :: x,y,z,vx,vy,vz
        integer::i
        ! remember, we are looping over the files. Therefore the same time is written to each file
        do i=1,N
            ! write on the file
            write(FILEUNITBASE+i,*)SNGL(myt),SNGL(x(i)),SNGL(y(i)),SNGL(z(i)),SNGL(vx(i)),SNGL(vy(i)),SNGL(vz(i))
        END DO


    END SUBROUTINE writeparticleorbits

    SUBROUTINE leapfrogintime(nstep,NP,xt,yt,zt,vxt,vyt,vzt)
        INTEGER, intent(in) :: nstep,NP ! number of time steps
        REAL*8, DIMENSION(NP,nstep+1), INTENT(OUT) :: xt,yt,zt,vxt,vyt,vzt

        REAL*8, DIMENSION(NP) :: ax,ay,az
        ! get the intermediate positions and velocities
        REAL*8, DIMENSION(NP) :: xtmp, ytmp, ztmp
        
        REAL*8 :: TESCTHRESHOLD = -999.0
        INTEGER :: i 
        INTEGER, DIMENSION(NP) :: indexes
        LOGICAL, DIMENSION(NP) :: isescaper
        ! for finding the energy with respect to the host and updating the escape time
        REAL*8, DIMENSION(NP) :: vx2host,vy2host,vz2host,Energy

        call assert_gravitational_constant_initialized()
        ! give each particle an index
        do i = 1,NP
            indexes(i) = i
        end do
        ! reset the index 
        i=0
        ! initalize the accelerations at zero
        ax = 0.0
        ay = 0.0
        az = 0.0
        ! initialize the positions and velocities
        xt=0
        yt=0
        zt=0
        vxt=0
        vyt=0
        vzt=0

        xt(:,1) = xf
        yt(:,1) = yf
        zt(:,1) = zf
        vxt(:,1) = vxf
        vyt(:,1) = vyf
        vzt(:,1) = vzf

        ! hit everything once to initialize the accelerations
        ! this is also necessary to initialize the host index
        currenttime = timestamps(1)
        call HIT(nparticles,xf,yf,zf,ax,ay,az)

        ! check if anyone is unbound 
        if (DOHOSTPERTURBER) then
            ! measure the energy of the particles with respect to the host
            vx2host = vxt(:,1)-host_vx_current
            vy2host = vyt(:,1)-host_vy_current
            vz2host = vzt(:,1)-host_vz_current
            Energy = 0.5*(vx2host**2+vy2host**2+vz2host**2) + phiHP
            ! update the escape time
            isescaper=(tesc < TESCTHRESHOLD .and. Energy> 0.0)
            tesc(PACK(indexes,isescaper)) = currenttime
        end if

        
        DO i= 1,nstep 
            currenttime = timestamps(i)
            ! drift a half step 
            xtmp = xt(:,i) + 0.5*dt*vxt(:, i)
            ytmp = yt(:,i) + 0.5*dt*vyt(:, i)
            ztmp = zt(:,i) + 0.5*dt*vzt(:, i)
            ! compute the accelerations at the initial time
            currenttime = (timestamps(i+1) + timestamps(i)) / 2.0
            call HIT(NP,xtmp,ytmp,ztmp,ax,ay,az)
            ! update the velocities a full step 
            vxt(:,i+1) = vxt(:,i) + ax*dt
            vyt(:,i+1) = vyt(:,i) + ay*dt
            vzt(:,i+1) = vzt(:,i) + az*dt   
            ! drift a half step 
            xt(:,i+1) = xtmp + 0.5*dt*vxt(:,i+1)
            yt(:,i+1) = ytmp + 0.5*dt*vyt(:,i+1)
            zt(:,i+1) = ztmp + 0.5*dt*vzt(:,i+1)
            currenttime = timestamps(i+1)


            if (DOHOSTPERTURBER) then
                vx2host = vxt(:,i+1)-host_vx_current
                vy2host = vyt(:,i+1)-host_vy_current
                vz2host = vzt(:,i+1)-host_vz_current
                Energy = 0.5*(vx2host**2+vy2host**2+vz2host**2) + phiHP
                ! update the escape time
                isescaper=(tesc < TESCTHRESHOLD .and. Energy> 0.0)
                tesc(PACK(indexes,isescaper)) = currenttime    
            end if
            
        END DO 
        ! store the final positions and velocities
        xf = xt(:,nstep+1)
        yf = yt(:,nstep+1)
        zf = zt(:,nstep+1)
        vxf = vxt(:,nstep+1)
        vyf = vyt(:,nstep+1)
        vzf = vzt(:,nstep+1)


    END SUBROUTINE  leapfrogintime

    SUBROUTINE leapfrogtofinalpositions()
        ! take the current positions and integrate until the end
        REAL*8, DIMENSION(nparticles) :: ax,ay,az
        ! REAL*8, DIMENSION(nparticles) :: xtemp,ytemp,ztemp 
        REAL*8 :: TESCTHRESHOLD = -999.0
        INTEGER :: i
        integer, dimension(nparticles) :: indexes
        logical, dimension(nparticles) :: isescaper
        ! for finding the energy with repsect to the host and updating the escape time
        REAL*8, DIMENSION(nparticles) :: vx2host,vy2host,vz2host,Energy 
        
        call assert_gravitational_constant_initialized()

        ! give each particle an index
        do i = 1,nparticles
            indexes(i) = i
        end do

        ! hit everything once to initialize the accelerations
        ! this is also necessary to initialize the host index
        currenttime = timestamps(1)
        call HIT(nparticles,xf,yf,zf,ax,ay,az)
        
        ! evaluate the potential at the initial positions
        if (DOHOSTPERTURBER) then
            ! measure the energy of the particles with respect to the host
            vx2host = vxf-host_vx_current
            vy2host = vyf-host_vy_current
            vz2host = vzf-host_vz_current
            Energy = 0.5*(vx2host**2+vy2host**2+vz2host**2) + phiHP
            ! update the escape time
            isescaper=(tesc < TESCTHRESHOLD .and. Energy> 0.0)
            tesc(PACK(indexes,isescaper)) = currenttime
        end if

        IF (DOWRITEORBITS) then
            CALL writeparticleorbits(currenttime,nparticles,xf,yf,zf,vxf,vyf,vzf)
        END IF
        if (DOWRITESTREAM) then
            CALL writestream(0,nparticles,xf,yf,zf,vxf,vyf,vzf)
        end if
        

        DO i=1,ntimesteps 
            currenttime = timestamps(i)
            ! first half drift 
            xf = xf + 0.5 * dt * vxf
            yf = yf + 0.5 * dt * vyf
            zf = zf + 0.5 * dt * vzf
            currenttime = (timestamps(i+1) + timestamps(i)) / 2.0
            ! compute the accelerations at the initial time
            call HIT(nparticles,xf,yf,zf,ax,ay,az)
            ! update the velocities a full step 
            vxf = vxf + ax*dt
            vyf = vyf + ay*dt
            vzf = vzf + az*dt   
            ! drift a half step 
            xf = xf +  0.5 * dt * vxf
            yf = yf +  0.5 * dt * vyf
            zf = zf +  0.5 * dt * vzf
            currenttime = timestamps(i+1)


            if (DOHOSTPERTURBER) then
                vx2host = vxf-host_vx_current
                vy2host = vyf-host_vy_current
                vz2host = vzf-host_vz_current
                Energy = 0.5d0*(vx2host**2+vy2host**2+vz2host**2) + phiHP
                ! update the escape time
                isescaper=(tesc < TESCTHRESHOLD .and. Energy> 0.0)
                tesc(PACK(indexes,isescaper)) = currenttime
            end if

            IF (DOWRITEORBITS) then
                if (MOD(i,nwriteskip).eq.0) then 
                    CALL writeparticleorbits(currenttime,nparticles,xf,yf,zf,vxf,vyf,vzf)
                end if 
            END IF    
            if (DOWRITESTREAM) then
                if (MOD(i,nwriteskip).eq.0) then 
                    CALL writestream(i/nwriteskip,nparticles,xf,yf,zf,vxf,vyf,vzf)
                end if 
            end if        
        END DO


    END SUBROUTINE leapfrogtofinalpositions


    SUBROUTINE ruthforestintime(nstep,NP,xt,yt,zt,vxt,vyt,vzt)
        ! integrate the positions and velocities forward in time
        ! return the positions and velocities at each timestep to the user
        INTEGER, intent(in) :: nstep,NP ! number of time steps
        REAL*8, DIMENSION(NP,nstep+1), INTENT(OUT) :: xt,yt,zt,vxt,vyt,vzt
        ! initialize the accelerations
        REAL*8, DIMENSION(NP) :: axf,ayf,azf
        REAL*8 :: TESCTHRESHOLD = -999.0
        INTEGER :: i
        integer, dimension(NP) :: indexes
        logical, dimension(NP) :: isescaper
        ! for finding the energy with repsect to the host and updating the escape time
        REAL*8, DIMENSION(NP) :: vx2host,vy2host,vz2host,Energy 
        ! THE coefficients for the Ruth-Forest integrator Waltz
        REAL*8  :: c1,c2,c3,c4,d1,d2,d3,d4 ! c for the positions, d for the velocities
        REAL*8  :: w ! for convience for coefficients
        INTEGER :: integration_sign

        call assert_gravitational_constant_initialized()

        if (DOBACKWARDORBIT) then
            integration_sign = -1
        else
            integration_sign = 1
        end if

        w = sqrt(2.0D0**(1.0D0/3.0D0) + 2.0D0**(-1.0D0/3.0D0) -1.0D0 )/6.0D0 ! D0 is for double precision

        c1 =  w + 0.5D0
        c2 = -w
        c3 = -w
        c4 =  w + 0.5D0

        d1 =  2.0D0*w+1.0D0
        d2 = -4.0D0*w-1.0D0
        d3 =  2.0D0*w+1.0D0
        d4 =  0.0D0
        ! give each particle an index
        do i = 1,NP
            indexes(i) = i
        end do
        ! reset the index 
        i=0
        ! initalize the accelerations at zero
        axf = 0.0
        ayf = 0.0
        azf = 0.0
        ! initialize the positions and velocities
        xt=0
        yt=0
        zt=0
        vxt=0
        vyt=0
        vzt=0
        xt(:,1) = xf
        yt(:,1) = yf
        zt(:,1) = zf
        vxt(:,1) = vxf
        vyt(:,1) = vyf
        vzt(:,1) = vzf

        currenttime=timestamps(1)
        ! hit everything once to initialize the accelerations
        ! this is also necessary to initialize the host index
        call HIT(nparticles,xf,yf,zf,axf,ayf,azf)
        
        if (DOHOSTPERTURBER) then
            ! measure the energy of the particles with respect to the host
            vx2host = vxf-host_vx_current
            vy2host = vyf-host_vy_current
            vz2host = vzf-host_vz_current
            Energy = 0.5*(vx2host**2+vy2host**2+vz2host**2) + phiHP
            ! update the escape time
            isescaper=(tesc < TESCTHRESHOLD .and. Energy> 0.0)
            tesc(PACK(indexes,isescaper)) = currenttime
        end if


        if (DEBUGACCELERATIONS) then
            call HIT(NP,xf,yf,zf,axf,ayf,azf)
            aSG(1,1) = axSG(1)
            aSG(2,1) = aySG(1)
            aSG(3,1) = azSG(1)
            aHP(1,1) = axHP(1)
            aHP(2,1) = ayHP(1)
            aHP(3,1) = azHP(1)
            aP(1,1) = axP(1)
            aP(2,1) = ayP(1)
            aP(3,1) = azP(1)
            aNBODY(1,1) = axNBODY(1)
            aNBODY(2,1) = ayNBODY(1)
            aNBODY(3,1) = azNBODY(1)
            aBAR(1,1) = axBAR(1)
            aBAR(2,1) = ayBAR(1)
            aBAR(3,1) = azBAR(1)
            aTOTAL(1,1) = axSG(1)+axHP(1)+axP(1)+axNBODY(1)+axBAR(1)
            aTOTAL(2,1) = aySG(1)+ayHP(1)+ayP(1)+ayNBODY(1)+ayBAR(1)
            aTOTAL(3,1) = azSG(1)+azHP(1)+azP(1)+azNBODY(1)+azBAR(1)
        end if 

        IF (DEBUGBARORIENTATION) then
            bartheta(1) = bar_angle
        end if

        do i=1,nstep
            currenttime=timestamps(i+1)
            ! drift
            xf = xf + c1*vxf*dt
            yf = yf + c1*vyf*dt
            zf = zf + c1*vzf*dt
            currenttime = currenttime + integration_sign*c1*dt
            ! kick
            call HIT(NP,xf,yf,zf,axf,ayf,azf)
            vxf = vxf + d1*axf*dt
            vyf = vyf + d1*ayf*dt
            vzf = vzf + d1*azf*dt
            ! drift
            xf = xf + c2*vxf*dt
            yf = yf + c2*vyf*dt
            zf = zf + c2*vzf*dt
            ! kick
            currenttime = currenttime + integration_sign*c2*dt
            call HIT(NP,xf,yf,zf,axf,ayf,azf)
            vxf = vxf + d2*axf*dt
            vyf = vyf + d2*ayf*dt
            vzf = vzf + d2*azf*dt
            ! drift
            xf = xf + c3*vxf*dt
            yf = yf + c3*vyf*dt
            zf = zf + c3*vzf*dt
            ! kick
            currenttime = currenttime + integration_sign*c3*dt
            call HIT(NP,xf,yf,zf,axf,ayf,azf)
            vxf = vxf + d3*axf*dt
            vyf = vyf + d3*ayf*dt
            vzf = vzf + d3*azf*dt
            ! drift
            xf = xf + c4*vxf*dt
            yf = yf + c4*vyf*dt
            zf = zf + c4*vzf*dt
            ! kick
            currenttime = currenttime + integration_sign*c4*dt
            call HIT(NP,xf,yf,zf,axf,ayf,azf)
            vxf = vxf + d4*axf*dt
            vyf = vyf + d4*ayf*dt
            vzf = vzf + d4*azf*dt
            xt(:,i+1) = xf
            yt(:,i+1) = yf
            zt(:,i+1) = zf
            vxt(:,i+1) = vxf
            vyt(:,i+1) = vyf
            vzt(:,i+1) = vzf

            if (DOHOSTPERTURBER) then
                ! measure the energy of the particles with respect to the host
                vx2host = vxf-host_vx_current
                vy2host = vyf-host_vy_current
                vz2host = vzf-host_vz_current
                Energy = 0.5*(vx2host**2+vy2host**2+vz2host**2) + phiHP
                ! update the escape time
                isescaper=(tesc < TESCTHRESHOLD .and. Energy> 0.0)
                tesc(PACK(indexes,isescaper)) = currenttime
            end if

            if (DEBUGACCELERATIONS) then
                aSG(1,i+1) = axSG(1)
                aSG(2,i+1) = aySG(1)
                aSG(3,i+1) = azSG(1)
                aHP(1,i+1) = axHP(1)
                aHP(2,i+1) = ayHP(1)
                aHP(3,i+1) = azHP(1)
                aP(1,i+1) = axP(1)
                aP(2,i+1) = ayP(1)
                aP(3,i+1) = azP(1)
                aNBODY(1,i+1) = axNBODY(1)
                aNBODY(2,i+1) = ayNBODY(1)
                aNBODY(3,i+1) = azNBODY(1)
                aBAR(1,i+1) = axBAR(1)
                aBAR(2,i+1) = ayBAR(1)
                aBAR(3,i+1) = azBAR(1)
                aTOTAL(1,i+1) = axSG(1)+axHP(1)+axP(1)+axNBODY(1)+axBAR(1)
                aTOTAL(2,i+1) = aySG(1)+ayHP(1)+ayP(1)+ayNBODY(1)+ayBAR(1)
                aTOTAL(3,i+1) = azSG(1)+azHP(1)+azP(1)+azNBODY(1)+azBAR(1)
            end if 

            IF (DEBUGBARORIENTATION) then
                bartheta(i+1) = bar_angle
            end if
        end do

    end subroutine ruthforestintime

    SUBROUTINE HIT(NP,x,y,z,ax,ay,az)
        ! Compute the net force on NP particles. Force-only: does not evaluate
        ! static gravity potential. Host/perturber/bar potentials are populated
        ! as module-level side effects (needed for escape-energy bookkeeping).
        INTEGER, INTENT(IN) :: NP
        REAL*8, DIMENSION(NP), INTENT(IN) :: x,y,z
        REAL*8, DIMENSION(NP), INTENT(OUT) :: ax,ay,az

        ! reset the accelerations to zero
        axSG = 0.0
        aySG = 0.0
        azSG = 0.0
        axHP = 0.0
        ayHP = 0.0
        azHP = 0.0
        axP = 0.0
        ayP = 0.0
        azP = 0.0
        axNBODY = 0.0
        ayNBODY = 0.0
        azNBODY = 0.0
        axBAR=0.0
        ayBAR=0.0
        azBAR=0.0
        phiHP=0.0
        phiP=0.0
        phiBAR=0.0

        if (grav_gravity_finalized) then
            call grav_evaluategravityforces(NP,x,y,z,axSG,aySG,azSG)
        end if

        if (DOHOSTPERTURBER) then
            CALL host_find_time_idx(currenttime)
            call host_compute_force(nparticles,x,y,z,axHP,ayHP,azHP,phiHP)
        end if

        if (DOPERTURBERS) then
            call pert_find_time_idx(currenttime)
            call pert_compute_force(nparticles,x,y,z,axP,ayP,azP,phiP)
        end if

        if (DONBODY) then
            call pot_nbody_plummers(nbodyparams,nparticles,x,y,z,axNBODY,ayNBODY,azNBODY,phiTensor)
        end if

        if (DOGALACTICBAR) then
            CALL bar_update(currenttime)
            call bar_force_eval(nparticles,x,y,z,axBAR,ayBAR,azBAR,phiBAR)
        end if

        if (DEBUGACCELERATIONS.eqv..TRUE.) then 
            print*, "SG: "
            print*, axSG(1),aySG(1),azSG(1)
            if (DOHOSTPERTURBER) then
                print*, "HP: "
                print*, axHP(1),ayHP(1),azHP(1)
            end if
            if (DOPERTURBERS) then
                print*, "P: "
                print*, axP(1),ayP(1),azP(1)        
            end if 
            if (DONBODY) then
                print*, "NBODY: "
                print*, axNBODY(1),ayNBODY(1),azNBODY(1)
            end if
            if (DOGALACTICBAR) then
                print*, "BAR: "
                print*, axBAR(1),ayBAR(1),azBAR(1)
            end if 
        end if

        ax=axSG+axHP+axP+axNBODY+axBAR
        ay=aySG+ayHP+ayP+ayNBODY+ayBAR
        az=azSG+azHP+azP+azNBODY+azBAR

    END SUBROUTINE HIT

    SUBROUTINE DEALLOCATE
        ! deallocate the arrays
        integer::i
        if (INITIALKINEMATICSSET) then 
            DEALLOCATE(xf,yf,zf,vxf,vyf,vzf,tesc)
            DEALLOCATE(axSG,aySG,azSG)
            DEALLOCATE(axHP,ayHP,azHP)
            DEALLOCATE(axP,ayP,azP)
            DEALLOCATE(axNBODY,ayNBODY,azNBODY)
            DEALLOCATE(axBAR,ayBAR,azBAR)
            DEALLOCATE(phiSG,phiHP,phiP,phiBAR,phiNBODY)
            DEALLOCATE(phiTensor)
            INITIALKINEMATICSSET = .FALSE.
        end if
        
        if (INTEGRATIONPARAMETERSSET) then
            deallocate(timestamps)
            INTEGRATIONPARAMETERSSET = .FALSE.
        end if 
        CALL cleargravitycomponents()
        
        if (DOHOSTPERTURBER) then
            CALL host_deallocate_mod
            DOHOSTPERTURBER=.FALSE.
        end if
        
        IF (DOPERTURBERS) then
            CALL pert_deallocate
            DOPERTURBERS=.FALSE.
        end if

        if (DOGALACTICBAR) THEN
            CALL bar_deallocate
            DOGALACTICBAR=.FALSE.
        end if 

        if (DONBODY) then
            DEALLOCATE(nbodyparams)
            DONBODY=.FALSE.
        end if

        IF (DOWRITEORBITS) then
            DOWRITEORBITS=.FALSE.
            DO i=1,nparticles
                close(FILEUNITBASE+i)
            END DO
        END IF
        
        IF (DOBACKWARDORBIT) then
            DOBACKWARDORBIT=.FALSE.
        END IF

        IF (DOWRITESTREAM) THEN
            DOWRITESTREAM=.FALSE.
        END IF

        if (DEBUGACCELERATIONS) then
            DEALLOCATE(aSG,aHP,aP,aNBODY,aBAR,aTOTAL)
            DEBUGACCELERATIONS=.FALSE.
        end if

        if (DEBUGBARORIENTATION) then 
            DEALLOCATE(bartheta)
            DEBUGBARORIENTATION=.FALSE.
        end if
    END SUBROUTINE DEALLOCATE

    SUBROUTINE finalizesimulator()
        ! Commit step for simulator configuration.
        ! Idempotently finalizes the gravity module if it has components but
        ! is not yet finalized. Hard error if no gravity components are registered.
        IF (.NOT. grav_gravity_finalized) THEN
            CALL grav_finalizegravity()
        END IF
    END SUBROUTINE finalizesimulator

    SUBROUTINE evaluatepotentials(NP, x, y, z, phi)
        ! Evaluate the total gravitational potential from all active physics modules.
        ! Intended for snapshot-time output, diagnostic checks, and escape-energy
        ! bookkeeping. Not called during integration substeps.
        ! NOTE: host/perturber/bar potential evaluation recomputes forces as a
        ! side effect (those modules compute force+potential together).
        INTEGER, INTENT(IN) :: NP
        REAL*8, DIMENSION(NP), INTENT(IN) :: x, y, z
        REAL*8, DIMENSION(NP), INTENT(OUT) :: phi
        INTEGER :: i, j

        phiSG = 0.0
        phiHP = 0.0
        phiP  = 0.0
        phiBAR = 0.0
        phiNBODY = 0.0

        if (grav_gravity_finalized) then
            call grav_evaluategravitypotential(NP, x, y, z, phiSG)
        end if

        if (DOHOSTPERTURBER) then
            CALL host_find_time_idx(currenttime)
            call host_compute_force(NP, x, y, z, axHP, ayHP, azHP, phiHP)
        end if

        if (DOPERTURBERS) then
            call pert_find_time_idx(currenttime)
            call pert_compute_force(NP, x, y, z, axP, ayP, azP, phiP)
        end if

        if (DONBODY) then
            call pot_nbody_plummers(nbodyparams, NP, x, y, z, axNBODY, ayNBODY, azNBODY, phiTensor)
            DO i = 1, NP
                DO j = i, NP
                    phiNBODY(i) = phiNBODY(i) + phiTensor(i, j)
                END DO
            END DO
        end if

        if (DOGALACTICBAR) then
            CALL bar_update(currenttime)
            call bar_force_eval(NP, x, y, z, axBAR, ayBAR, azBAR, phiBAR)
        end if

        phi = phiSG + phiHP + phiP + phiBAR + phiNBODY
    END SUBROUTINE evaluatepotentials


END MODULE simulator


