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
            print*, "ERROR: setgalacticbar must be called before setdebugbarorientation"
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
        INTEGER, intent(in) :: nstep,NP ! number of time steps
        ! return the positions and velocities at each timestep to the user
        REAL*8, DIMENSION(NP,nstep+1), INTENT(OUT) :: xt,yt,zt,vxt,vyt,vzt
        ! initialize the accelerations
        REAL*8, DIMENSION(NP) :: axSG,aySG,azSG ! static galaxy
        REAL*8, DIMENSION(NP) :: axHP,ayHP,azHP ! host perturber
        REAL*8, DIMENSION(NP) :: axP,ayP,azP ! perturbers
        REAL*8, DIMENSION(NP) :: axNBODY,ayNBODY,azNBODY ! nbody
        REAL*8, DIMENSION(NP) :: axBAR,ayBAR,azBAR ! bar
        REAL*8, DIMENSION(NP) :: axf,ayf,azf ! total
        REAL*8, DIMENSION(NP) :: ax0,ay0,az0
        REAL*8, DIMENSION(NP) :: phiSG,phiHP,phiP,phiBAR
        REAL*8, DIMENSION(NP,NP) :: phiNBODY
        REAL*8 :: TESCTHRESHOLD = -999.0
        INTEGER :: i,j
        integer, dimension(NP) :: indexes
        logical, dimension(NP) :: isescaper
        ! for finding the energy with repsect to the host and updating the escape time
        REAL*8, DIMENSION(NP) :: vx2host,vy2host,vz2host,Energy 
        ! give each particle an index
        do i = 1,NP
            indexes(i) = i
        end do
        ! initalize the accelerations at zero
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
        ax0 = 0.0
        ay0 = 0.0
        az0 = 0.0
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
        i=0
        xt(:,1) = xf
        yt(:,1) = yf
        zt(:,1) = zf
        vxt(:,1) = vxf
        vyt(:,1) = vyf
        vzt(:,1) = vzf
        ! compute the accelerations at the initial time
        call milkywaypotential(milkwayparams,NP,xt(:,1),yt(:,1),zt(:,1),axSG,aySG,azSG,phiSG)
        IF (DOGALACTICBAR) THEN
            CALL updatebarorientation(currenttime)
            CALL barforce(NP,xt(:,1),yt(:,1),zt(:,1),axBAR,ayBAR,azBAR,phiBAR)
            if (DEBUGBARORIENTATION) then
                bartheta(1) = theta
                timestamps(1)=currenttime
            end if            

        END IF
        if (DOHOSTPERTURBER) then
            CALL findhosttimeindex(currenttime)
            CALL computeforcebyhosts(NP,xt(:,1),yt(:,1),zt(:,1),axHP,ayHP,azHP,phiHP)
            ! measure the energy of the particles with respect to the host
            vx2host = vxt(:,1)-vxhost(hosttimeindex)
            vy2host = vyt(:,1)-vyhost(hosttimeindex)
            vz2host = vzt(:,1)-vzhost(hosttimeindex)
            Energy = 0.5*(vx2host**2+vy2host**2+vz2host**2) + phiHP
            ! update the escape time
            isescaper=(tesc < TESCTHRESHOLD .and. Energy> 0.0)
            tesc(PACK(indexes,isescaper)) = currenttime
        end if
        IF (DONBODY) then              
            CALL NBODYPLUMMERS(nbodyparams,NP,xt(:,1),yt(:,1),zt(:,1),axNBODY,ayNBODY,azNBODY,phiNBODY)
        end if 

        ax0=axSG+axHP+axP+axNBODY+axBAR
        ay0=aySG+ayHP+ayP+ayNBODY+ayBAR
        az0=azSG+azHP+azP+azNBODY+azBAR
        ! for debugging 
        if (DEBUGACCELERATIONS) then
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
            aTOTAL(1,1) = ax0(1)
            aTOTAL(2,1) = ay0(1)
            aTOTAL(3,1) = az0(1)
        end if
        DO i=1,(nstep)
            currenttime=timestamps(i)
            xt(:,i+1) = xt(:,i) + vxt(:,i)*dt + 0.5*ax0*dt**2
            yt(:,i+1) = yt(:,i) + vyt(:,i)*dt + 0.5*ay0*dt**2
            zt(:,i+1) = zt(:,i) + vzt(:,i)*dt + 0.5*az0*dt**2
            call milkywaypotential(milkwayparams,NP,xt(:,i+1),yt(:,i+1),zt(:,i+1),axSG,aySG,azSG,phiSG)
            if (DOHOSTPERTURBER) then
                CALL advancehosttimeindex()
                CALL computeforcebyhosts(NP,xt(:,i+1),yt(:,i+1),zt(:,i+1),axHP,ayHP,azHP,phiHP)
            end if 
            IF (DOGALACTICBAR) THEN
                CALL updatebarorientation(currenttime)
                CALL barforce(NP,xt(:,i+1),yt(:,i+1),zt(:,i+1),axBAR,ayBAR,azBAR,phiBAR)
                if (DEBUGBARORIENTATION) then
                    bartheta(i+1) = theta
                end if
            END IF
            IF (DONBODY) then
                CALL NBODYPLUMMERS(nbodyparams,NP,xt(:,i+1),yt(:,i+1),zt(:,i+1),axNBODY,ayNBODY,azNBODY,phiNBODY)
            end if     
            ! measure the energy of the particles with respect to the host
            axf=axSG+axHP+axP+axNBODY+axBAR
            ayf=aySG+ayHP+ayP+ayNBODY+ayBAR
            azf=azSG+azHP+azP+azNBODY+azBAR
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
                aTOTAL(1,i+1) = axf(1)
                aTOTAL(2,i+1) = ayf(1)
                aTOTAL(3,i+1) = azf(1)
            end if

        END DO
    END SUBROUTINE leapfrogintime

    SUBROUTINE leapfrogtofinalpositions()
        ! take the current positions and integrate until the end
        REAL*8, DIMENSION(nparticles) :: axSG,aySG,azSG ! static galaxy
        REAL*8, DIMENSION(nparticles) :: axHP,ayHP,azHP ! host perturber
        REAL*8, DIMENSION(nparticles) :: axP,ayP,azP ! perturbers
        REAL*8, DIMENSION(nparticles) :: axNBODY,ayNBODY,azNBODY ! nbody
        REAL*8, DIMENSION(nparticles) :: axBAR,ayBAR,azBAR ! bar
        REAL*8, DIMENSION(nparticles) :: axf,ayf,azf ! total
        REAL*8, DIMENSION(nparticles) :: ax0,ay0,az0
        REAL*8, DIMENSION(nparticles) :: phiSG,phiHP,phiP,phiBAR
        REAL*8, DIMENSION(nparticles,nparticles) :: phiNBODY
        REAL*8 :: TESCTHRESHOLD = -999.0
        INTEGER :: i,j
        integer, dimension(nparticles) :: indexes
        logical, dimension(nparticles) :: isescaper
        REAL*8, DIMENSION(nparticles) :: x0,y0,z0,vx0,vy0,vz0
        ! for finding the energy with repsect to the host and updating the escape time
        REAL*8, DIMENSION(nparticles) :: vx2host,vy2host,vz2host,Energy 
        
        print*, "leapfrogtofinalpositions"
        ! give each particle an index
        do i = 1,nparticles
            indexes(i) = i
        end do
        !! initalize the accelerations at zero
        axSG = 0.0
        aySG = 0.0
        azSG = 0.0
        axHP = 0.0
        ayHP = 0.0
        azHP = 0.0
        axP = 0.0
        ayP = 0.0
        azP = 0.0
        axBAR = 0.0
        ayBAR = 0.0
        azBAR = 0.0
        axNBODY = 0.0
        ayNBODY = 0.0
        azNBODY = 0.0
        ax0 = 0.0
        ay0 = 0.0
        az0 = 0.0
        axf = 0.0
        ayf = 0.0
        azf = 0.0
        ! set up the initial positions
        x0 = xf
        y0 = yf
        z0 = zf
        vx0 = vxf
        vy0 = vyf
        vz0 = vzf
        ! evaluate the potential at the initial positions
        call milkywaypotential(milkwayparams,nparticles,x0,y0,z0,axSG,aySG,azSG,phiSG)
        IF (DOGALACTICBAR) THEN
            CALL updatebarorientation(currenttime)
            CALL barforce(nparticles,x0,y0,z0,axBAR,ayBAR,azBAR,phiBAR)
        END IF
        if (DOPERTURBERS) THEN
            call findperturbertimeindex(currenttime)
            call computeforcebyperturbers(nparticles,x0,y0,z0,axP,ayP,azP,phiP)
        end IF
        if (DOHOSTPERTURBER) then
            CALL findhosttimeindex(currenttime)
            CALL computeforcebyhosts(nparticles,x0,y0,z0,axHP,ayHP,azHP,phiHP)
            ! measure the energy of the particles with respect to the host
            vx2host = vx0-vxhost(hosttimeindex)
            vy2host = vy0-vyhost(hosttimeindex)
            vz2host = vz0-vzhost(hosttimeindex)
            Energy = 0.5*(vx2host**2+vy2host**2+vz2host**2) + phiHP
            ! update the escape time
            isescaper=(tesc < TESCTHRESHOLD .and. Energy> 0.0)
            tesc(PACK(indexes,isescaper)) = currenttime
        end if

        IF (DONBODY) then
            CALL NBODYPLUMMERS(nbodyparams,nparticles,x0,y0,z0,axNBODY,ayNBODY,azNBODY,phiNBODY)
        end if     

        ax0=axSG+axHP+axP+axNBODY+axBAR
        ay0=aySG+ayHP+ayP+ayNBODY+ayBAR
        az0=azSG+azHP+azP+azNBODY+azBAR
        IF (DOWRITEORBITS) then
            CALL writeparticleorbits(currenttime,nparticles,x0,y0,z0,vx0,vy0,vz0)
        END IF
        if (DOWRITESTREAM) then
            CALL writestream(0,nparticles,x0,y0,z0,vx0,vy0,vz0)
        end if
        DO i=1,(ntimesteps)
            currenttime=timestamps(i)
            xf = x0 + vx0*dt + 0.5*ax0*dt**2
            yf = y0 + vy0*dt + 0.5*ay0*dt**2
            zf = z0 + vz0*dt + 0.5*az0*dt**2
            call milkywaypotential(milkwayparams,nparticles,xf,yf,zf,axSG,aySG,azSG,phiSG)
            IF (DOGALACTICBAR) THEN
                CALL updatebarorientation(currenttime)
                CALL barforce(nparticles,x0,y0,z0,axBAR,ayBAR,azBAR,phiBAR)
            END IF    
            if (DOHOSTPERTURBER) then
                CALL advancehosttimeindex()
                CALL computeforcebyhosts(nparticles,xf,yf,zf,axHP,ayHP,azHP,phiHP)
            end if 
            if (DOPERTURBERS) THEN
                call advanceperturbertimeindex(currenttime)
                call computeforcebyperturbers(nparticles,x0,y0,z0,axP,ayP,azP,phiP)
            end IF            
            IF (DONBODY) then
                CALL NBODYPLUMMERS(nbodyparams,nparticles,xf,yf,zf,axNBODY,ayNBODY,azNBODY,phiNBODY)
            end if     
            ! measure the energy of the particles with respect to the host
            axf=axSG+axHP+axP+axNBODY+axBAR
            ayf=aySG+ayHP+ayP+ayNBODY+ayBAR
            azf=azSG+azHP+azP+azNBODY+axBAR
            vxf = vx0 + 0.5*(ax0+axf)*dt
            vyf = vy0 + 0.5*(ay0+ayf)*dt
            vzf = vz0 + 0.5*(az0+azf)*dt   
            ! update the positions and velocities
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


    SUBROUTINE DEALLOCATE
        ! deallocate the arrays
        integer::i
        if (INITIALKINEMATICSSET) then 
            DEALLOCATE(xf,yf,zf,vxf,vyf,vzf,tesc)
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


