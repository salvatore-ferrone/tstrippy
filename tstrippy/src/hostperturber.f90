MODULE hostperturber
    ! so FAR, we can either have constant mass, or a double exponential mass evolution
    
    use potentials, only : plummer
    use mathutils, only : linear_interp_scalar

    IMPLICIT NONE

    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: xhost,yhost,zhost,vxhost,vyhost,vzhost,timehost
    REAL*8, PUBLIC :: masshost=0.0D0, radiushost=0.0D0
    REAL*8, PUBLIC :: xhostcurrent=0.0D0,yhostcurrent=0.0D0,zhostcurrent=0.0D0
    REAL*8, PUBLIC :: vxhostcurrent=0.0D0,vyhostcurrent=0.0D0,vzhostcurrent=0.0D0
    REAL*8, PUBLIC :: masshostcurrent=0.0D0, radiushostcurrent=0.0D0
    ! timehost: must be an ordered list from smallest to largest (negative to positive)
    INTEGER, PUBLIC :: hosttimeindex = 1
    PUBLIC :: host_init_kinematics, host_init_mass, host_init_radius
    PUBLIC :: findhosttimeindex, updatehoststate
    PUBLIC :: hostallocation, hostdeallocation, computeforcebyhosts

    INTEGER, PUBLIC :: host_mass_model = 0 ! 0=constant, 1=double exponential, 2=polynomial, etc.
    REAL*8, PUBLIC, ALLOCATABLE :: host_mass_params(:)

    CONTAINS
    
    SUBROUTINE host_init_kinematics(NTIMESTEPS, t, x, y, z, vx, vy, vz)
        INTEGER, INTENT(IN) :: NTIMESTEPS
        REAL*8, INTENT(IN), DIMENSION(NTIMESTEPS) :: t, x, y, z, vx, vy, vz
        INTEGER :: i

        if (NTIMESTEPS < 2) then
            print*, "ERROR: host kinematics need at least 2 time points for interpolation"
            stop
        end if

        do i=2,NTIMESTEPS
            if (t(i) <= t(i-1)) then
                print*, "ERROR: host times must be strictly increasing"
                stop
            end if
        end do

        CALL hostallocation(NTIMESTEPS)
        xhost = x
        yhost = y
        zhost = z
        vxhost = vx
        vyhost = vy
        vzhost = vz
        timehost = t
        hosttimeindex = 1
    END SUBROUTINE host_init_kinematics

    SUBROUTINE host_init_mass(mass_model, params)
        INTEGER, INTENT(IN) :: mass_model
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        host_mass_model = mass_model
        IF (ALLOCATED(host_mass_params)) DEALLOCATE(host_mass_params)
        IF (host_mass_model == 0) THEN
            masshost = params(1)
            ALLOCATE(host_mass_params(1))
            host_mass_params(1) = params(1)
        ELSE
            ALLOCATE(host_mass_params(SIZE(params)))
            host_mass_params = params
        END IF
    END SUBROUTINE host_init_mass

    SUBROUTINE host_init_radius(radius)
        REAL*8, INTENT(IN) :: radius
        radiushost = radius
        radiushostcurrent = radius
    END SUBROUTINE host_init_radius

    FUNCTION get_host_mass(t)
        REAL*8 :: t, get_host_mass
        SELECT CASE (host_mass_model)
        CASE (0)
            get_host_mass = host_mass_params(1) ! constant mass
        CASE (1)
            get_host_mass = mass_double_exponential(t, host_mass_params)
        ! Add more CASEs for new models
        CASE DEFAULT
            get_host_mass = host_mass_params(1) ! fallback
        END SELECT
    END FUNCTION get_host_mass

    FUNCTION mass_double_exponential(t, params)
        REAL*8 :: t, params(:), mass_double_exponential
        REAL*8 :: A1, tau1, A2, tau2, Mf, tf, C
        A1 = params(1)
        tau1 = params(2)
        A2 = params(3)
        tau2 = params(4)
        Mf = params(5)
        tf = params(6)
        C = Mf - (A1 * EXP(-tf/tau1) + A2 * EXP(-tf/tau2))
        mass_double_exponential = A1 * EXP(-t/tau1) + A2 * EXP(-t/tau2) + C
    END FUNCTION mass_double_exponential    

    SUBROUTINE updatehoststate(mytime)
        REAL*8, INTENT(IN) :: mytime
        REAL*8 :: t0, t1, alpha, denom
        INTEGER :: ntime

        ntime = SIZE(timehost)
        if (ntime < 2) then
            print*, "ERROR: host time array has fewer than 2 points"
            stop
        end if

        if (mytime <= timehost(1)) then
            hosttimeindex = 1
            alpha = 0.0D0
        else if (mytime >= timehost(ntime)) then
            hosttimeindex = ntime - 1
            alpha = 1.0D0
        else
            ! Track the lower bracketing index around mytime.
            DO WHILE (hosttimeindex < ntime - 1 .AND. timehost(hosttimeindex + 1) < mytime)
                hosttimeindex = hosttimeindex + 1
            END DO

            DO WHILE (hosttimeindex > 1 .AND. timehost(hosttimeindex) > mytime)
                hosttimeindex = hosttimeindex - 1
            END DO

            t0 = timehost(hosttimeindex)
            t1 = timehost(hosttimeindex + 1)
            denom = t1 - t0
            if (denom <= 0.0D0) then
                print*, "ERROR: host times are not strictly increasing"
                stop
            end if
            alpha = (mytime - t0) / denom
        end if

        xhostcurrent = linear_interp_scalar(xhost(hosttimeindex), xhost(hosttimeindex + 1), alpha)
        yhostcurrent = linear_interp_scalar(yhost(hosttimeindex), yhost(hosttimeindex + 1), alpha)
        zhostcurrent = linear_interp_scalar(zhost(hosttimeindex), zhost(hosttimeindex + 1), alpha)
        vxhostcurrent = linear_interp_scalar(vxhost(hosttimeindex), vxhost(hosttimeindex + 1), alpha)
        vyhostcurrent = linear_interp_scalar(vyhost(hosttimeindex), vyhost(hosttimeindex + 1), alpha)
        vzhostcurrent = linear_interp_scalar(vzhost(hosttimeindex), vzhost(hosttimeindex + 1), alpha)

        masshost = get_host_mass(mytime)
        masshostcurrent = masshost
        radiushostcurrent = radiushost
    END SUBROUTINE updatehoststate

    SUBROUTINE findhosttimeindex(mytime)
        REAL*8, INTENT(IN) :: mytime
        ! Backward-compatible entry point used by the integrator.
        CALL updatehoststate(mytime)
    END SUBROUTINE findhosttimeindex

    SUBROUTINE computeforcebyhosts(Nparticles,Gin,x,y,z,ax,ay,az,phi)
        ! compute the force on the particles due to the hosts
        ! the force is computed by summing over all hosts
        ! the force is computed in galactic coordinates
        real*8, intent(in) :: Gin
        integer, intent(in) :: Nparticles
        real*8, intent(in), dimension(Nparticles) :: x,y,z
        real*8, intent(out), dimension(Nparticles) :: ax,ay,az,phi
        REAL*8, dimension(3) :: params
        real*8,dimension(Nparticles) :: dx,dy,dz,axhost,ayhost,azhost,phihost
        params(1) =  Gin
        ax = 0
        ay = 0
        az = 0
        phi = 0

        dx = x - xhostcurrent
        dy = y - yhostcurrent
        dz = z - zhostcurrent
        params(2) = masshostcurrent
        params(3) = radiushostcurrent
        call plummer(params,Nparticles,dx,dy,dz,axhost,ayhost,azhost,phihost)
        ax = ax + axhost
        ay = ay + ayhost
        az = az + azhost
        phi = phi + phihost
    END SUBROUTINE computeforcebyhosts

    subroutine hostallocation(NTIMESTEPS)
        integer, intent(in) ::  NTIMESTEPS
        allocate(xhost(NTIMESTEPS))
        allocate(yhost(NTIMESTEPS))
        allocate(zhost(NTIMESTEPS))
        allocate(vxhost(NTIMESTEPS))
        allocate(vyhost(NTIMESTEPS))
        allocate(vzhost(NTIMESTEPS))
        allocate(timehost(NTIMESTEPS))
    END SUBROUTINE hostallocation

    SUBROUTINE hostdeallocation
        if (allocated(xhost)) deallocate(xhost)
        if (allocated(yhost)) deallocate(yhost)
        if (allocated(zhost)) deallocate(zhost)
        if (allocated(vxhost)) deallocate(vxhost)
        if (allocated(vyhost)) deallocate(vyhost)
        if (allocated(vzhost)) deallocate(vzhost)
        if (allocated(timehost)) deallocate(timehost)
        if (allocated(host_mass_params)) deallocate(host_mass_params)
    END SUBROUTINE hostdeallocation    

END MODULE hostperturber



