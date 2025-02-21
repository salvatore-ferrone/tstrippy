MODULE perturbers
    ! This module contains the perturbers of the system
    ! For now, the perturbers are other plummers
    ! all perturbers have the same timesteps
    ! so we have NPERTURBERS * Ntimesteps for x,y,z 
    ! then we have a one dimensional array of the time
    ! everything should be in galactic coordinates
    use potentials, only : plummer
    IMPLICIT NONE
    PRIVATE
    REAL*8, DIMENSION(:,:), PUBLIC, ALLOCATABLE :: xperturbers,yperturbers,zperturbers 
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: timeperturbers,massperturber,radiusperturber
    ! timeperturbers: must be an ordered list from smallest to largest (negative to positive)
    INTEGER, PUBLIC :: perturbertimeindex 
    PUBLIC :: perturberinitialization,findperturbertimeindex,advanceperturbertimeindex
    PUBLIC :: perturberallocation,perturberdeallocation,computeforcebyperturbers
    REAL*8,parameter :: G=4.300917270036279e-06 !! in solar masses and km/s
    CONTAINS
    
    ! initialize the perturbers
    subroutine perturberinitialization(NPERTURBERS,NTIMESTEPS,t,x,y,z,mass,radius)
        integer, intent(in) :: NPERTURBERS, NTIMESTEPS
        real*8, intent(in), dimension(NPERTURBERS) :: mass,radius
        real*8, intent(in), dimension(NPERTURBERS,NTIMESTEPS) :: x,y,z
        real*8, intent(in), dimension(NTIMESTEPS):: t

        call perturberallocation(NPERTURBERS,NTIMESTEPS)
        xperturbers = x
        yperturbers = y
        zperturbers = z
        timeperturbers = t
        massperturber = mass
        radiusperturber =radius

        
    end subroutine perturberinitialization

    subroutine perturberallocation(NPERTURBERS,NTIMESTEPS)
        integer, intent(in) :: NPERTURBERS, NTIMESTEPS
        allocate(xperturbers(NPERTURBERS,NTIMESTEPS))
        allocate(yperturbers(NPERTURBERS,NTIMESTEPS))
        allocate(zperturbers(NPERTURBERS,NTIMESTEPS))
        allocate(timeperturbers(NTIMESTEPS))
        allocate(massperturber(NPERTURBERS))
        allocate(radiusperturber(NPERTURBERS))
    END SUBROUTINE perturberallocation


    SUBROUTINE findperturbertimeindex(mytime)
        ! find the index of the perturber that is just below mytime
        real*8, intent(in) :: mytime
        real*8:: dt,globalmin
        INTEGER :: i
        globalmin=abs(mytime-timeperturbers(1))
        DO i=1,size(timeperturbers)
            dt=abs(mytime-timeperturbers(i))
            IF (dt <= globalmin) THEN
                perturbertimeindex = i
            END IF
        END DO
    END SUBROUTINE findperturbertimeindex

    SUBROUTINE advanceperturbertimeindex(mytime)
        ! make sure the time index of the system is just above mytime. 
        ! if it is not, advance it until it is
        ! also, set an upper limit such that the time index is never larger than the number of timesteps
        real*8, intent(in) :: mytime
        DO while ((mytime > timeperturbers(perturbertimeindex)))
            perturbertimeindex = perturbertimeindex + 1
            if (perturbertimeindex > size(timeperturbers)) then
                perturbertimeindex = size(timeperturbers)
                exit
            end if
        END DO 
    END SUBROUTINE advanceperturbertimeindex



    SUBROUTINE computeforcebyperturbers(Nparticles,x,y,z,ax,ay,az,phi)
        ! compute the force on the particles due to the perturbers
        ! the force is computed by summing over all perturbers
        ! the force is computed in galactic coordinates
        integer, intent(in) :: Nparticles
        real*8, intent(in), dimension(Nparticles) :: x,y,z
        real*8, intent(out), dimension(Nparticles) :: ax,ay,az,phi
        REAL*8, dimension(3) :: params
        real*8,dimension(Nparticles) :: dx,dy,dz,axperturber,ayperturber,azperturber,phiperturber
        integer :: i,nperturbers
        nperturbers=size(massperturber)
        params(1) = G
        ax = 0
        ay = 0
        az = 0
        phi = 0
        DO i=1,nperturbers
            dx = x - xperturbers(i,perturbertimeindex)
            dy = y - yperturbers(i,perturbertimeindex)
            dz = z - zperturbers(i,perturbertimeindex)
            params(2) = massperturber(i)
            params(3) = radiusperturber(i)
            call plummer(params,Nparticles,dx,dy,dz,axperturber,ayperturber,azperturber,phiperturber)
            ax = ax + axperturber
            ay = ay + ayperturber
            az = az + azperturber
            phi = phi + phiperturber
        END DO 
    END SUBROUTINE computeforcebyperturbers
    ! deallocate the perturbers
    SUBROUTINE perturberdeallocation
        deallocate(xperturbers,yperturbers,zperturbers,timeperturbers,massperturber,radiusperturber)
    END SUBROUTINE perturberdeallocation

END MODULE perturbers



