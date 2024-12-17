MODULE integrator 
    ! this integrator needs to contain the current positions
    ! it needs to be able to apply any force that I want at any time
    ! it needs to be able to integrate the positions and velocities
    use constants, only : G
    use potentials
    use perturbers
    use hostperturber
    use galacticbar
    IMPLICIT NONE
    PRIVATE 
    ! DECLARE SUBROUTINES
    PUBLIC :: setstaticgalaxy,setintegrationparameters,setinitialkinematics
    PUBLIC :: setdebugaccelerations, setdebugbarorientation,setbackwardorbit
    PUBLIC :: inithostperturber,initnbodysystem,initgalacticbar,initperturbers
    PUBLIC :: leapfrogtofinalpositions,leapfrogintime
    PUBLIC :: ruthforestintime
    PUBLIC :: HIT
    PUBLIC :: initwriteparticleorbits,writeparticleorbits
    PUBLIC :: initwritestream,writestream
    PUBLIC :: deallocate
    ! DECIDE WHICH PHYSICS TO INCLUDE
    LOGICAL, PUBLIC :: DONBODY = .FALSE.
    LOGICAL, PUBLIC :: DOPERTURBERS = .FALSE.
    LOGICAL, PUBLIC :: DOHOSTPERTURBER = .FALSE.
    LOGICAL, PUBLIC :: DOGALACTICBAR = .FALSE.
    LOGICAL, PUBLIC :: DOBACKWARDORBIT = .FALSE.
    ! Variables to keep track of the physics that has been set
    LOGICAL, PUBLIC :: GALAXYISSET = .FALSE.
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
    procedure(), pointer,public :: milkywaypotential
    REAL*8,DIMENSION(:),PUBLIC,allocatable :: milkwayparams
    REAL*8, PUBLIC :: currenttime,dt
    INTEGER, PUBLIC :: ntimesteps,ntimepoints,nparticles,nwriteskip
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
    SUBROUTINE setstaticgalaxy(milkywaypotentialname,mwparams)
        !! decide which potential we are going to integrate in
        !! more error checking will be done at the python level
        character*100, intent(in) :: milkywaypotentialname
        REAL*8, DIMENSION(:), intent(in) :: mwparams
        INTEGER :: nparams
        if (milkywaypotentialname.EQ."pouliasis2017pii") then
            milkywaypotential => pouliasis2017PII
        else if (milkywaypotentialname.EQ."miyamotonagai") then
            milkywaypotential => miyamotonagai
        else if (milkywaypotentialname.EQ."allensantillianhalo") then
            milkywaypotential => allensantillianhalo
        else if (milkywaypotentialname.EQ."plummer") then
            milkywaypotential => plummer
        else if (milkywaypotentialname.EQ."longmuralibar") then
            milkywaypotential => longmuralibar
        else
            print*, "ERROR. milkywaypotential not found"
            print*, "the string must be a valid potential name from potentials.f90"
            stop
        end if
        nparams = size(mwparams)
        allocate(milkwayparams(nparams))
        milkwayparams = mwparams
        GALAXYISSET=.TRUE.
    END SUBROUTINE setstaticgalaxy


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
    
    subroutine initnbodysystem(N,massesnbody,scaleradiinbody)
        ! initialize the nbody system
        ! meaning that everytime the system is evaluated, we also compute the Nbody forces
        INTEGER, intent(in) :: N
        REAL*8, DIMENSION(N), intent(in) :: massesnbody,scaleradiinbody
        DONBODY = .TRUE.
        allocate(nbodyparams(2*N+1))
        nbodyparams(1)=G
        nbodyparams(1+1:N+1)=massesnbody
        nbodyparams(N+1+1:2*N+1)=scaleradiinbody
    end subroutine initnbodysystem

    subroutine inithostperturber(nhosttimepoints,timeH,xH,yH,zH,vxH,vyH,vzH,massh,radiush)
        ! initialize the host perturber
        INTEGER, intent(in) :: nhosttimepoints
        real*8, intent(in) :: massh,radiush
        real*8, intent(in), dimension(nhosttimepoints) :: timeH,xH,yH,zH,vxH,vyH,vzH
        DOHOSTPERTURBER = .TRUE.
        CALL hostinitialization(nhosttimepoints,timeH,xH,yH,zH,vxH,vyH,vzH,massh,radiush)
    end subroutine inithostperturber

    SUBROUTINE initperturbers(nperturbers,nperturbertimesteps,tp,xp,yp,zp,masses,radii)
        integer, intent(in) :: nperturbers,nperturbertimesteps
        real*8, intent(in) ,dimension(nperturbertimesteps) :: tp
        real*8, intent(in), dimension(nperturbers,nperturbertimesteps) :: xp,yp,zp
        real*8, intent(in), dimension(nperturbers) :: masses,radii
        DOPERTURBERS = .TRUE.
        CALL perturberinitialization(nperturbers,nperturbertimesteps,tp,xp,yp,zp,masses,radii)

    END SUBROUTINE initperturbers

    SUBROUTINE initgalacticbar(barpotenname,barparams,barpoly)
        ! initialize the galactic bar
        character*100, INTENT(IN) :: barpotenname
        REAL*8, DIMENSION(:), intent(in) :: barparams,barpoly
        DOGALACTICBAR = .TRUE.
        CALL galacticbarinitialization(barpotenname,barparams,barpoly)
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
        ! integrate the positions and velocities forward in time
        ! return the positions and velocities at each timestep to the user
        INTEGER, intent(in) :: nstep,NP ! number of time steps
        REAL*8, DIMENSION(NP,nstep+1), INTENT(OUT) :: xt,yt,zt,vxt,vyt,vzt
        ! initialize the accelerations
        REAL*8, DIMENSION(NP) :: axf,ayf,azf 
        REAL*8, DIMENSION(NP) :: ax0,ay0,az0
        REAL*8, DIMENSION(NP) :: phi
        REAL*8 :: TESCTHRESHOLD = -999.0
        INTEGER :: i
        integer, dimension(NP) :: indexes
        logical, dimension(NP) :: isescaper
        ! for finding the energy with repsect to the host and updating the escape time
        REAL*8, DIMENSION(NP) :: vx2host,vy2host,vz2host,Energy 
        ! give each particle an index
        do i = 1,NP
            indexes(i) = i
        end do
        ! reset the index 
        i=0
        ! initalize the accelerations at zero
        ax0 = 0.0
        ay0 = 0.0
        az0 = 0.0
        axf = 0.0
        ayf = 0.0
        azf = 0.0
        phi = 0.0
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
        ! compute the accelerations at the initial time
        call HIT(NP,xt(:,1),yt(:,1),zt(:,1),ax0,ay0,az0,phi)
        ! check for unbound particles
        if (DOHOSTPERTURBER) then
            ! measure the energy of the particles with respect to the host
            vx2host = vxt(:,1)-vxhost(hosttimeindex)
            vy2host = vyt(:,1)-vyhost(hosttimeindex)
            vz2host = vzt(:,1)-vzhost(hosttimeindex)
            Energy = 0.5*(vx2host**2+vy2host**2+vz2host**2) + phiHP
            ! update the escape time
            isescaper=(tesc < TESCTHRESHOLD .and. Energy> 0.0)
            tesc(PACK(indexes,isescaper)) = currenttime
        end if


        DO i=1,(nstep)
            currenttime=timestamps(i+1)
            xt(:,i+1) = xt(:,i) + vxt(:,i)*dt + 0.5*ax0*dt**2
            yt(:,i+1) = yt(:,i) + vyt(:,i)*dt + 0.5*ay0*dt**2
            zt(:,i+1) = zt(:,i) + vzt(:,i)*dt + 0.5*az0*dt**2
            call HIT(NP,xt(:,i+1),yt(:,i+1),zt(:,i+1),axf,ayf,azf,phi)
            vxt(:,i+1) = vxt(:,i) + 0.5*(ax0+axf)*dt
            vyt(:,i+1) = vyt(:,i) + 0.5*(ay0+ayf)*dt
            vzt(:,i+1) = vzt(:,i) + 0.5*(az0+azf)*dt   
            ax0=axf
            ay0=ayf
            az0=azf
            if (DOHOSTPERTURBER) then
                vx2host = vxt(:,i+1)-vxhost(hosttimeindex)
                vy2host = vyt(:,i+1)-vyhost(hosttimeindex)
                vz2host = vzt(:,i+1)-vzhost(hosttimeindex)
                Energy = 0.5*(vx2host**2+vy2host**2+vz2host**2) + phiHP
                ! update the escape time
                isescaper=(tesc < TESCTHRESHOLD .and. Energy> 0.0)
                tesc(PACK(indexes,isescaper)) = currenttime    
            end if


        END DO
    END SUBROUTINE leapfrogintime

    SUBROUTINE leapfrogtofinalpositions()
        ! take the current positions and integrate until the end
        REAL*8, DIMENSION(nparticles) :: axf,ayf,azf,phi ! total
        REAL*8, DIMENSION(nparticles) :: ax0,ay0,az0

        REAL*8 :: TESCTHRESHOLD = -999.0
        INTEGER :: i
        integer, dimension(nparticles) :: indexes
        logical, dimension(nparticles) :: isescaper
        REAL*8, DIMENSION(nparticles) :: x0,y0,z0,vx0,vy0,vz0
        ! for finding the energy with repsect to the host and updating the escape time
        REAL*8, DIMENSION(nparticles) :: vx2host,vy2host,vz2host,Energy 
        

        ! give each particle an index
        do i = 1,nparticles
            indexes(i) = i
        end do
        ! set up the initial positions
        x0 = xf
        y0 = yf
        z0 = zf
        vx0 = vxf
        vy0 = vyf
        vz0 = vzf
        ! evaluate the potential at the initial positions


        call HIT(nparticles,x0,y0,z0,ax0,ay0,az0,phi)
        
        if (DOHOSTPERTURBER) then
            ! measure the energy of the particles with respect to the host
            vx2host = vx0-vxhost(hosttimeindex)
            vy2host = vy0-vyhost(hosttimeindex)
            vz2host = vz0-vzhost(hosttimeindex)
            Energy = 0.5*(vx2host**2+vy2host**2+vz2host**2) + phiHP
            ! update the escape time
            isescaper=(tesc < TESCTHRESHOLD .and. Energy> 0.0)
            tesc(PACK(indexes,isescaper)) = currenttime
        end if

        IF (DOWRITEORBITS) then
            CALL writeparticleorbits(currenttime,nparticles,x0,y0,z0,vx0,vy0,vz0)
        END IF
        if (DOWRITESTREAM) then
            CALL writestream(0,nparticles,x0,y0,z0,vx0,vy0,vz0)
        end if
        
        DO i=1,ntimesteps
            currenttime=timestamps(i+1)
            xf = x0 + vx0*dt + 0.5*ax0*dt**2
            yf = y0 + vy0*dt + 0.5*ay0*dt**2
            zf = z0 + vz0*dt + 0.5*az0*dt**2
            call HIT(nparticles,xf,yf,zf,axf,ayf,azf,phi)
            vxf = vx0 + 0.5*(ax0+axf)*dt
            vyf = vy0 + 0.5*(ay0+ayf)*dt
            vzf = vz0 + 0.5*(az0+azf)*dt   

            x0=xf
            y0=yf
            z0=zf
            vx0=vxf
            vy0=vyf
            vz0=vzf
            ax0=axf
            ay0=ayf
            az0=azf
            if (DOHOSTPERTURBER) then
                vx2host = vx0-vxhost(hosttimeindex)
                vy2host = vy0-vyhost(hosttimeindex)
                vz2host = vz0-vzhost(hosttimeindex)
                Energy = 0.5*(vx2host**2+vy2host**2+vz2host**2) + phiHP
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
        REAL*8, DIMENSION(NP) :: phi
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
        phi = 0.0
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
        if (DOHOSTPERTURBER) then
            ! measure the energy of the particles with respect to the host
            vx2host = vxf-vxhost(hosttimeindex)
            vy2host = vyf-vyhost(hosttimeindex)
            vz2host = vzf-vzhost(hosttimeindex)
            Energy = 0.5*(vx2host**2+vy2host**2+vz2host**2) + phiHP
            ! update the escape time
            isescaper=(tesc < TESCTHRESHOLD .and. Energy> 0.0)
            tesc(PACK(indexes,isescaper)) = currenttime
        end if


        if (DEBUGACCELERATIONS) then
            call HIT(NP,xf,yf,zf,axf,ayf,azf,phi)
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
            bartheta(1) = theta
        end if
        do i=1,nstep
            currenttime=timestamps(i)
            ! drift
            xf = xf + c1*vxf*dt
            yf = yf + c1*vyf*dt
            zf = zf + c1*vzf*dt
            currenttime = currenttime + integration_sign*c1*dt
            ! kick
            call HIT(NP,xf,yf,zf,axf,ayf,azf,phi)
            vxf = vxf + d1*axf*dt
            vyf = vyf + d1*ayf*dt
            vzf = vzf + d1*azf*dt
            ! drift
            xf = xf + c2*vxf*dt
            yf = yf + c2*vyf*dt
            zf = zf + c2*vzf*dt
            ! kick
            currenttime = currenttime + integration_sign*c2*dt
            call HIT(NP,xf,yf,zf,axf,ayf,azf,phi)
            vxf = vxf + d2*axf*dt
            vyf = vyf + d2*ayf*dt
            vzf = vzf + d2*azf*dt
            ! drift
            xf = xf + c3*vxf*dt
            yf = yf + c3*vyf*dt
            zf = zf + c3*vzf*dt
            ! kick
            currenttime = currenttime + integration_sign*c3*dt
            call HIT(NP,xf,yf,zf,axf,ayf,azf,phi)
            vxf = vxf + d3*axf*dt
            vyf = vyf + d3*ayf*dt
            vzf = vzf + d3*azf*dt
            ! drift
            xf = xf + c4*vxf*dt
            yf = yf + c4*vyf*dt
            zf = zf + c4*vzf*dt
            ! kick
            currenttime = currenttime + integration_sign*c4*dt
            call HIT(NP,xf,yf,zf,axf,ayf,azf,phi)
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
                vx2host = vxf-vxhost(hosttimeindex)
                vy2host = vyf-vyhost(hosttimeindex)
                vz2host = vzf-vzhost(hosttimeindex)
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
                bartheta(i+1) = theta
            end if
        end do

    end subroutine ruthforestintime

    SUBROUTINE HIT(NP,x,y,z,ax,ay,az,phi)
        ! compute the acceleration of the particles at a given position 
        INTEGER, INTENT(IN) :: NP
        REAL*8, DIMENSION(NP), INTENT(IN) :: x,y,z
        REAL*8, DIMENSION(NP), INTENT(OUT) :: ax,ay,az
        REAL*8, DIMENSION(NP), INTENT(OUT) :: phi
        INTEGER :: i,j ! loop variables for summing over phiTensor

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
        phiSG=0.0
        phiHP=0.0
        phiP=0.0
        phiBAR=0.0
        phiNBODY=0.0
        
        if (GALAXYISSET) then
            call milkywaypotential(milkwayparams,nparticles,x,y,z,axSG,aySG,azSG,phiSG)
        end if

        if (DOHOSTPERTURBER) then
            CALL findhosttimeindex(currenttime)
            call computeforcebyhosts(nparticles,x,y,z,axHP,ayHP,azHP,phiHP)
        end if

        if (DOPERTURBERS) then
            call findperturbertimeindex(currenttime)
            call computeforcebyperturbers(nparticles,x,y,z,axP,ayP,azP,phiP)
        end if

        if (DONBODY) then
            call NBODYPLUMMERS(nbodyparams,nparticles,x,y,z,axNBODY,ayNBODY,azNBODY,phiTensor)
        end if

        if (DOGALACTICBAR) then
            CALL updatebarorientation(currenttime)
            call barforce(nparticles,x,y,z,axBAR,ayBAR,azBAR,phiBAR)
        end if

        ax=axSG+axHP+axP+axNBODY+axBAR
        ay=aySG+ayHP+ayP+ayNBODY+ayBAR
        az=azSG+azHP+azP+azNBODY+azBAR

        DO i=1,NP
            DO j=i,NP
                phiNBODY(i) = phiNBODY(i) + phiTensor(i,j)
            END DO
        END DO
        
        phi=phiSG+phiHP+phiP+phiBAR+phiNBODY
        

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
        IF (GALAXYISSET) then
            DEALLOCATE(milkwayparams)
            GALAXYISSET = .FALSE.
        END IF 
        
        if (DOHOSTPERTURBER) then
            CALL hostdeallocation
            DOHOSTPERTURBER=.FALSE.
        end if
        
        IF (DOPERTURBERS) then
            CALL perturberdeallocation
            DOPERTURBERS=.FALSE.
        end if

        if (DOGALACTICBAR) THEN
            CALL bardeallocation
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


END MODULE integrator


