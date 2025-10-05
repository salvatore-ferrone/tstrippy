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
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: timeperturbers
    REAL*8, DIMENSION(:,:),PUBLIC, ALLOCATABLE:: massperturber, radiusperturber
    ! timeperturbers: must be an ordered list from smallest to largest (negative to positive)
    INTEGER, PUBLIC :: perturbertimeindex 
    PUBLIC :: perturberinitialization,findperturbertimeindex,advanceperturbertimeindex
    PUBLIC :: perturberdeallocation,computeforcebyperturbers
    PRIVATE :: taylor_eval
    REAL*8, PUBLIC :: G
    CONTAINS
    
    ! initialize the perturbers
    subroutine perturberinitialization(t,x,y,z,Gin,mass,radius)
        real*8, intent(in), dimension(:,:) :: mass,radius
        real*8, intent(in), dimension(:,:) :: x,y,z
        real*8, intent(in), dimension(:):: t
        real*8, intent(in) :: Gin
        INTEGER :: NPERTURBERS, NTIMESTEPS, NMASSCOEFF, NRADIUSCOEFF
        
        ! check that we are the same size
        if (size(x,1) /= size(y,1) .or. size(x,1) /= size(z,1)) then
            print *, "Error: x, y, z must have the same first dimension size!"
            stop
        end if
        ! same with the mass and radius
        if (size(x,1) /= size(mass,1)) then
            print *, "Error: mass must have the same first dimension size as x, y, z!"
            stop
        end if
        ! if the radius is not the same size as the mass, stop
        if (size(x,1) /= size(radius,1)) then
            print *, "Error: radius must have the same first dimension size as x, y, z!"
            stop
        end if
        ! Make sure the number of timesteps are good too 
        if (size(x,2) /= size(y,2) .or. size(x,2) /= size(z,2)) then
            print *, "Error: x, y, z must have the same second dimension size!"
            stop
        end if
        
        if (size(t) /= size(x,2)) then
            print *, "Error: t must have the same size as the second dimension of x, y, z!"
            stop
        end if         
        
        ! store the mass and radius coefficents
        NPERTURBERS = size(x,1)
        NTIMESTEPS = size(x,2)
        NMASSCOEFF = size(mass,2)
        NRADIUSCOEFF = size(radius,2)
        

        ! allocate them 
        allocate(xperturbers(NPERTURBERS,NTIMESTEPS))
        allocate(yperturbers(NPERTURBERS,NTIMESTEPS))
        allocate(zperturbers(NPERTURBERS,NTIMESTEPS))
        allocate(timeperturbers(NTIMESTEPS))
        ALLOCATE(massperturber(NPERTURBERS,NMASSCOEFF))
        ALLOCATE(radiusperturber(NPERTURBERS,NRADIUSCOEFF))

        xperturbers = x
        yperturbers = y
        zperturbers = z
        timeperturbers = t
        massperturber = mass
        radiusperturber = radius
        G = Gin
    end subroutine perturberinitialization




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

    FUNCTION taylor_eval(coefficients, t) RESULT(myvalue)
        REAL*8, INTENT(IN), DIMENSION(:) :: coefficients
        REAL*8, INTENT(IN) :: t
        REAL*8 :: myvalue
        integer :: ncoefficients, i
        ncoefficients = SIZE(coefficients)
        myvalue = 0.0D0
        DO i = 1, ncoefficients
            myvalue = myvalue + coefficients(i) * t**(i-1)
        END DO
    END FUNCTION taylor_eval

    SUBROUTINE computeforcebyperturbers(Nparticles,x,y,z,ax,ay,az,phi)
        ! compute the force on the particles due to the perturbers
        ! the force is computed by summing over all perturbers
        ! the force is computed in galactic coordinates
        integer, intent(in) :: Nparticles
        real*8, intent(in), dimension(Nparticles) :: x,y,z
        real*8, intent(out), dimension(Nparticles) :: ax,ay,az,phi
        REAL*8, dimension(3) :: params
        real*8,dimension(Nparticles) :: dx,dy,dz,axperturber,ayperturber,azperturber,phiperturber
        REAL*8 ::  current_time
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
            current_time = timeperturbers(perturbertimeindex)
            params(2) = taylor_eval(massperturber(i,:), current_time)
            params(3) = taylor_eval(radiusperturber(i,:), current_time)
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



