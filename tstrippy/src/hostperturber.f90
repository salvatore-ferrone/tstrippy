MODULE hostperturber
    ! so FAR, we can either have constant mass, or a double exponential mass evolution
    
    use potentials, only : plummer

    IMPLICIT NONE

    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: xhost,yhost,zhost,vxhost,vyhost,vzhost,timehost
    REAL*8, PUBLIC :: masshost, radiushost
    ! timehost: must be an ordered list from smallest to largest (negative to positive)
    INTEGER, PUBLIC :: hosttimeindex = 1
    PUBLIC :: host_init_kinematics, host_init_mass, host_init_radius
    PUBLIC :: findhosttimeindex
    PUBLIC :: hostallocation, hostdeallocation, computeforcebyhosts

    INTEGER, PUBLIC :: host_mass_model = 0 ! 0=constant, 1=double exponential, 2=polynomial, etc.
    REAL*8, PUBLIC, ALLOCATABLE :: host_mass_params(:)

    CONTAINS
    
    ! Kinematics initialization
    SUBROUTINE host_init_kinematics(NTIMESTEPS, t, x, y, z, vx, vy, vz)
        INTEGER, INTENT(IN) :: NTIMESTEPS
        REAL*8, INTENT(IN), DIMENSION(NTIMESTEPS) :: t, x, y, z, vx, vy, vz
        CALL hostallocation(NTIMESTEPS)
        xhost = x
        yhost = y
        zhost = z
        vxhost = vx
        vyhost = vy
        vzhost = vz
        timehost = t
    END SUBROUTINE host_init_kinematics

    ! Mass initialization (constant or model)
    SUBROUTINE host_init_mass(mass_model, params)
        INTEGER, INTENT(IN) :: mass_model
        REAL*8, INTENT(IN), DIMENSION(:) :: params
        host_mass_model = mass_model
        IF (host_mass_model == 0) THEN
            masshost = params(1)
            ALLOCATE(host_mass_params(1))
            host_mass_params(1) = params(1)
        ELSE
            IF (ALLOCATED(host_mass_params)) DEALLOCATE(host_mass_params)
            ALLOCATE(host_mass_params(SIZE(params)))
            host_mass_params = params
        END IF
    END SUBROUTINE host_init_mass

    ! Radius initialization
    SUBROUTINE host_init_radius(radius)
        REAL*8, INTENT(IN) :: radius
        radiushost = radius
    END SUBROUTINE host_init_radius




    ! Generic mass evolution function
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

    subroutine hostallocation(NTIMESTEPS)
        integer, intent(in) ::  NTIMESTEPS
        allocate(xhost(NTIMESTEPS))
        allocate(yhost(NTIMESTEPS))
        allocate(zhost(NTIMESTEPS))
        allocate(timehost(NTIMESTEPS))
    END SUBROUTINE hostallocation

    ! deallocate the hosts
    SUBROUTINE hostdeallocation
        deallocate(xhost,yhost,zhost,timehost)
    END SUBROUTINE hostdeallocation

    SUBROUTINE findhosttimeindex(mytime)
        ! This searche is done to find the closest time index in the timehost array
        ! It doesn't use MINLOC like it did before because that is an O(N) operation
        ! now we take advantage of the current index and move forward or backward
        ! until we find the closest time index to mytime
        ! This is an O(1) operation in the best case and O(N) in the worst case
        ! but it is much faster than the previous implementation
        ! Then we perform a check to see if the next index is closer
        ! this is important because if we are using leapfrog, we need to make sure
        ! that we are using the middle index to ensure time-reversability


        REAL*8, INTENT(IN) :: mytime
        INTEGER :: next_idx
        REAL*8 :: dist_current, dist_next
        
        ! Use the current hosttimeindex as starting point for search
        ! Move forward if needed
        DO WHILE (hosttimeindex < size(timehost) .AND. timehost(hosttimeindex) < mytime)
            hosttimeindex = hosttimeindex + 1
        END DO
        
        ! Move backward if we went too far
        DO WHILE (hosttimeindex > 1 .AND. timehost(hosttimeindex) > mytime)
            hosttimeindex = hosttimeindex - 1
        END DO
        
        ! Now find which is closer - current index or next index
        IF (hosttimeindex < size(timehost)) THEN
            next_idx = hosttimeindex + 1
            dist_current = ABS(timehost(hosttimeindex) - mytime)
            dist_next = ABS(timehost(next_idx) - mytime)
            
            ! Choose the closer timestamp
            IF (dist_next < dist_current) THEN
                hosttimeindex = next_idx
            END IF
        END IF
        
        ! update the host mass
        masshost = get_host_mass(timehost(hosttimeindex))
    
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

        dx = x - xhost(hosttimeindex)
        dy = y - yhost(hosttimeindex)
        dz = z - zhost(hosttimeindex)
        params(2) = masshost
        params(3) = radiushost
        call plummer(params,Nparticles,dx,dy,dz,axhost,ayhost,azhost,phihost)
        ax = ax + axhost
        ay = ay + ayhost
        az = az + azhost
        phi = phi + phihost

    END SUBROUTINE computeforcebyhosts

END MODULE hostperturber



