MODULE galacticbar
    ! this is the module for the galactic bar
    ! so far it will just include one mode
    ! it will handle updating the bar as a function of time
    ! it will return the force from the bar at a given time
		!    test comment from github app 
    ! The bar can have a changing angular momentum, 
    ! the orientation is found with a polynomial fit
    ! theta = theta0 + b*t + c*t^2 + d*t^3 + ...
    ! THUS  barorientationpolynomailcoeffs(1) = theta0
    ! use constants
    use potentials
    IMPLICIT NONE
    PRIVATE 
    ! DECLARE SUBROUTINES
    PUBLIC :: galacticbarinitialization,updatebarorientation,transformtobarframe
    PUBLIC :: transformtogalacticframe,barforce,bardeallocation
    ! DECLATE VARIABLES
    REAL*8, dimension(:), PUBLIC, ALLOCATABLE :: barpotentialparameters
    REAL*8, dimension(:), PUBLIC, ALLOCATABLE :: barorientationpolynomailcoeffs
    REAL*8, PUBLIC :: theta0,theta
    procedure(),pointer,public :: barpotential
    REAL*8, PARAMETER :: pi = 2.0d0 * acos(0.0d0)
    CONTAINS 
    SUBROUTINE galacticbarinitialization(barpotentialname,barpotenparams,barpolycoeff)
        ! this subroutine initializes the galactic bar
        character*100, intent(in) :: barpotentialname
        REAL*8, dimension(:), intent(in) :: barpotenparams
        REAL*8, dimension(:), intent(in) :: barpolycoeff
        INTEGER :: npotentialparameters, npolycoeffs
        IF (barpotentialname.EQ."longmuralibar") then 
            barpotential => longmuralibar
            npotentialparameters = 5
            ! make sure that the number of parameters is correct
            IF (size(barpotenparams).NE.npotentialparameters) then
                WRITE(*,*) "ERROR: INCORRECT NUMBER OF PARAMETERS FOR LONGMURALIBAR"
                STOP
            END IF
        ELSE 
            WRITE(*,*) "ERROR: INCORRECT BAR POTENTIAL NAME. must be a bar from potentials.f90"
            STOP
        END IF
        npolycoeffs=size(barpolycoeff)
        ! allocate the arrays
        ALLOCATE(barpotentialparameters(npotentialparameters))
        ALLOCATE(barorientationpolynomailcoeffs(npolycoeffs))
        ! set the values
        barpotentialparameters=barpotenparams
        barorientationpolynomailcoeffs=barpolycoeff
        theta0=barpolycoeff(1)
        theta=barpolycoeff(1)
    END SUBROUTINE galacticbarinitialization


    SUBROUTINE updatebarorientation(time)
        ! find the current orientation of the bar
        ! this is done with a polynomial fit of the orientation
        ! theta = theta0 + b*t + c*t^2 + d*t^3 + ...
        REAL*8, intent(in) :: time
        REAL*8 :: thetaTemp
        INTEGER :: i
        thetaTemp=0
        DO i=1,size(barorientationpolynomailcoeffs)
            thetaTemp=thetaTemp+barorientationpolynomailcoeffs(i)*time**(i-1)
        END DO
        theta=MOD(thetaTemp,2*pi)
        ! print*, "time, theta = ", time,theta
    END SUBROUTINE updatebarorientation



    SUBROUTINE transformtobarframe(N,xp,yp,zp,xout,yout,zout)
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN),DIMENSION(N):: xp,yp,zp
        REAL*8, INTENT(OUT),DIMENSION(N) :: xout,yout,zout
        xout=xp*cos(theta)+yp*sin(theta)
        yout=-xp*sin(theta)+yp*cos(theta)
        zout=zp
    END SUBROUTINE transformtobarframe

    SUBROUTINE transformtogalacticframe(N,xp,yp,zp,xout,yout,zout)
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN),DIMENSION(N):: xp,yp,zp
        REAL*8, INTENT(OUT),DIMENSION(N) :: xout,yout,zout
        xout=xp*cos(-theta)+yp*sin(-theta)
        yout=-xp*sin(-theta)+yp*cos(-theta)
        zout=zp    
    END SUBROUTINE transformtogalacticframe

    SUBROUTINE barforce(N,xp,yp,zp,ax,ay,az,phi)
        ! this subroutine calculates the force from the bar
        ! the input is the position of the particle in the galactic frame
        ! the output is the acceleration in the galactic frame
        ! the acceleration is calculated in the bar frame and then transformed to the galactic frame
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN),DIMENSION(N):: xp,yp,zp
        REAL*8, INTENT(OUT),DIMENSION(N) :: ax,ay,az
        REAL*8, INTENT(OUT),DIMENSION(N) :: phi
        REAL*8, DIMENSION(N) :: xbar,ybar,zbar
        REAL*8, DIMENSION(N) :: axbar,aybar,azbar
        CALL transformtobarframe(N,xp,yp,zp,xbar,ybar,zbar)
        CALL barpotential(barpotentialparameters,N,xbar,ybar,zbar,axbar,aybar,azbar,phi)
        CALL transformtogalacticframe(N,axbar,aybar,azbar,ax,ay,az) ! is it correct to transform the acceleration?
    END SUBROUTINE barforce

    SUBROUTINE bardeallocation()
        DEALLOCATE(barpotentialparameters)
        DEALLOCATE(barorientationpolynomailcoeffs)
    END SUBROUTINE bardeallocation
END MODULE galacticbar
