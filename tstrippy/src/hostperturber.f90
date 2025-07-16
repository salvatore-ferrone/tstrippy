MODULE hostperturber
    ! this is the exact same class as pertuerbers. 
    ! this is because I don't think F90 does inheritence 
    ! also I do not think F90 can do two instances of the same class in the same program
    use potentials, only : plummer

    IMPLICIT NONE

    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: xhost,yhost,zhost,vxhost,vyhost,vzhost,timehost
    REAL*8, PUBLIC :: G, masshost, radiushost
    ! timehost: must be an ordered list from smallest to largest (negative to positive)
    INTEGER, PUBLIC :: hosttimeindex = 1
    PUBLIC :: hostinitialization,findhosttimeindex,advancehosttimeindex
    PUBLIC :: hostallocation,hostdeallocation,computeforcebyhosts
    ! REAL*8,parameter :: G=4.300917270036279e-06 !! in solar masses and km/s
    CONTAINS
    
    ! initialize the hosts
    subroutine hostinitialization(NTIMESTEPS,t,x,y,z,vx,vy,vz,Gin,mass,radius)
        ! the radius is the plummer radius
        integer, intent(in) ::  NTIMESTEPS
        real*8, intent(in) :: mass,radius
        real*8, intent(in) :: Gin
        real*8, intent(in), dimension(NTIMESTEPS) :: x,y,z,vx,vy,vz
        real*8, intent(in), dimension(NTIMESTEPS):: t
        call hostallocation(NTIMESTEPS)
        xhost = x
        yhost = y
        zhost = z
        vxhost=vx
        vyhost=vy
        vzhost=vz
        timehost = t
        G = Gin
        masshost = mass
        radiushost =radius
        
    end subroutine hostinitialization

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
    END SUBROUTINE findhosttimeindex

    SUBROUTINE advancehosttimeindex()
        ! make sure the time index of the system is just above mytime. 
        ! if it is not, advance it until it is
        ! also, set an upper limit such that the time index is never larger than the number of timesteps
        hosttimeindex=hosttimeindex+1
        ! DO while ((mytime > timehost(hosttimeindex)))
        !     hosttimeindex = hosttimeindex + 1
        !     if (hosttimeindex > size(timehost)) then
        !         hosttimeindex = size(timehost)
        !         exit
        !     end if
        ! END DO 
    END SUBROUTINE advancehosttimeindex


    SUBROUTINE computeforcebyhosts(Nparticles,x,y,z,ax,ay,az,phi)
        ! compute the force on the particles due to the hosts
        ! the force is computed by summing over all hosts
        ! the force is computed in galactic coordinates
        integer, intent(in) :: Nparticles
        real*8, intent(in), dimension(Nparticles) :: x,y,z
        real*8, intent(out), dimension(Nparticles) :: ax,ay,az,phi
        REAL*8, dimension(3) :: params
        real*8,dimension(Nparticles) :: dx,dy,dz,axhost,ayhost,azhost,phihost
        params(1) = G
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



