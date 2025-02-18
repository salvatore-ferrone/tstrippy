!python -m numpy.f2py -c temp.f90 -m temp

! adding another comment!!
MODULE temp 
    ! USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE 

    PROCEDURE(), pointer, public :: my_system
    REAL*8, PARAMETER :: PI = 2.0d0 * acos(0.0d0)
    contains 



    FUNCTION KingDensityW(W) RESULT(rho)
        ! King density profile as a function of W
        ! W is the minimium of the potential divided by the velocity dispersion squared
        REAL*8, INTENT(IN) :: W 
        REAL*8 :: rho 
        REAL*8 :: sqW 
        sqW = SQRT(W)
        rho = EXP(W) * erf_custom(W) - (2.0/SQRT(PI))*sqW * (1.0 + 2.0*W/3.0)
    end FUNCTION

    FUNCTION KingScaleRadius(W0) result(r0)
        REAL*8, INTENT(IN) :: W0
        REAL*8 :: rho0
        REAL*8 :: r0
        rho0=KingDensityW(W0)
        r0 = sqrt(9.0 / (4.0 * PI * rho0) )
    end FUNCTION

    FUNCTION king_break_radius(W0) result(rbreak)
        ! used to separate the core and the halo when solving the ODE
        ! when solving, we will be linear in r when r< rbreak
        ! and we will sample w linearly when r > rbreak
        REAL*8, INTENT(IN) :: W0
        REAL*8 :: rbreak
        REAL*8 :: r0, threshold
        r0=KingScaleRadius(W0)
        threshold=2.0
        rbreak = r0 
        if (W0.lt.threshold) then
            rbreak = r0/100.0
        end if
    end FUNCTION

    SUBROUTINE king_ode_in_r(t,y,dydt)
        REAL*8, INTENT(IN) :: t
        REAL*8, INTENT(IN),dimension(2) :: y
        REAL*8, INTENT(OUT),dimension(2) :: dydt
        REAL*8 :: r,W, dwdr, d2wdr2, rho 
        ! unpack the variables 
        r = t
        W = y(1)
        dwdr=y(2)
        ! numerical gaurd 
        if (r.le.0.0) then 
            d2wdr2=0.0
        else 
            rho = KingDensityW(W)
            d2wdr2 = -4.0 * PI * rho - 2*dwdr/r
        end if 
        dydt(1) = dwdr
        dydt(2) = d2wdr2
    END SUBROUTINE king_ode_in_r

    SUBROUTINE king_ode_in_w(t,y,dydt)
        REAL*8, INTENT(IN) :: t
        REAL*8, INTENT(IN),dimension(2) :: y
        REAL*8, INTENT(OUT),dimension(2) :: dydt
        REAL*8 :: r,W, dwdr, d2wdr2, rho, term1, term2
        ! unpack the variables, in this case W is the independent variable 
        W = t
        r = y(1)
        dwdr=y(2)
        rho = KingDensityW(W)
        term1 = -1.0*(1.0/dwdr)
        term2 = (4.0 * PI * rho + 2.0*dwdr/r)
        print*, rho, term1, term2
        d2wdr2 = term1 * term2
        dydt(1) = 1.0/dwdr
        dydt(2) = d2wdr2
    END SUBROUTINE king_ode_in_w



    subroutine rk4(system_name,t_span,y0,nparams,params,npoints,nvars,tout,yout)
        REAL*8, intent(in), dimension(nvars):: y0
        REAL*8, intent(in), dimension(2):: t_span
        INTEGER, INTENT(in) :: npoints,nparams,nvars
        REAL*8, intent(out), dimension(npoints) :: tout
        REAL*8, intent(out), dimension(npoints,nvars) :: yout 
        REAL*8, dimension(nparams), intent(IN):: params
        REAL*8 :: dt,t0,tf
        REAL*8, dimension(nvars) :: k1,k2,k3,k4,y_temp
        character*100, intent(in) :: system_name

        ! coefficients for the coefficients 
        integer :: i
        t0=t_span(1)
        tf=t_span(2)
        if (system_name.eq."exponential_decay") then 
            my_system=>exponential_decay
        else if (system_name.eq."double_system_of_equations") then 
            my_system=>double_system_of_equations
        else
            print*, "error, system of equations not implemented"
            stop 
        end if 
        ! compute the time step 
        dt = (tf-t0) / (npoints - 1 )
        ! set the initial condition
        yout(1,:) = y0
        tout(1) = t0
        DO i = 2, npoints
            ! update the equation
            call my_system(tout(i-1),yout(i-1,:),k1,params)
            y_temp =  yout(i-1,:) + 0.5d0*dt*k1
            call my_system(tout(i-1) + 0.5d0*dt, y_temp, k2, params)
            y_temp = yout(i-1,:) + 0.5d0*dt*k2
            call my_system(tout(i-1) + 0.5d0*dt, y_temp, k3, params)
            y_temp = yout(i-1,:) + dt*k3
            call my_system(tout(i-1) + dt, y_temp, k4, params)
            yout(i,:) = yout(i-1,:) + (k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4) * dt / 6.0d0
            tout(i) = tout(i-1) + dt
        END DO 
    end subroutine rk4


        
    FUNCTION erf_custom(x) RESULT(y)
        ! Abramowitz and Stegun approximation for the error function provided to me by co-pilot
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: x
        REAL*8 :: y
        REAL*8, PARAMETER :: p = 0.3275911d0
        REAL*8, PARAMETER :: a1 = 0.254829592d0, a2 = -0.284496736d0
        REAL*8, PARAMETER :: a3 = 1.421413741d0, a4 = -1.453152027d0, a5 = 1.061405429d0
        REAL*8 :: t
        t = 1.0d0 / (1.0d0 + p * abs(x))
        y = 1.0d0 - ((( ( (a5 * t + a4) * t + a3) * t + a2) * t + a1) * t) * exp(-x * x)
        IF (x < 0.0d0) y = -y
    END FUNCTION erf_custom        


    ! dummy functions for building the module
    subroutine exponential_decay(t,y_in,dydt,params)
        REAL*8, INTENT(IN) :: t,y_in
        REAL*8, INTENT(OUT) :: dydt
        REAL*8, intent(in), dimension(1) :: params
        REAL*8 :: coeff
        coeff = params(1)
        dydt = -y_in*coeff
    end subroutine exponential_decay

    subroutine double_system_of_equations(t,y_in,dydt,params)
        REAL*8, INTENT(IN) :: t
        REAL*8, INTENT(IN),dimension(2) :: y_in
        REAL*8, INTENT(OUT),dimension(2) :: dydt
        REAL*8, intent(in), dimension(1) :: params
        REAL*8 :: coeff
        coeff = params(1)
        dydt(1) = y_in(2)
        dydt(2) = -y_in(1)*coeff
    end subroutine double_system_of_equations


end module temp 
    