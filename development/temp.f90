!python -m numpy.f2py -c temp.f90 -m temp

MODULE temp 
    IMPLICIT NONE 

    REAL*8 :: W0_
    INTEGER :: npoints_
    REAL*8, DIMENSION(:), ALLOCATABLE :: r_, W_, dwdr_
    PROCEDURE(), pointer, public :: my_system
    REAL*8, PARAMETER :: PI = 2.0d0 * acos(0.0d0)
    contains 

    SUBROUTINE initialize_king_potential_profile(W0, npoints)
        REAL*8, INTENT(IN) :: W0
        INTEGER, INTENT(IN) :: npoints
        REAL*8, ALLOCATABLE :: r(:), W(:), dwdr(:)

        ! Allocate local arrays
        ALLOCATE(r(npoints), W(npoints), dwdr(npoints))

        ! Call the solve_king_potential_profile subroutine
        CALL solve_king_potential_profile(W0, npoints, r, W, dwdr)

        ! Store the results in global variables
        W0_ = W0
        ALLOCATE(r_(npoints), W_(npoints), dwdr_(npoints))
        npoints_ = npoints
        r_ = r
        W_ = W
        dwdr_ = dwdr

        ! Deallocate local arrays
        DEALLOCATE(r, W, dwdr)
    END SUBROUTINE initialize_king_potential_profile

    SUBROUTINE solve_king_potential_profile(W0,npoints,r,W,dwdr)
        REAL*8, INTENT(IN) :: W0
        INTEGER, INTENT(IN) :: npoints
        REAL*8, INTENT(OUT), DIMENSION(npoints) :: r, W, dwdr
        REAL*8 :: rbreak 
        REAL*8, DIMENSION(2) :: y0,t_span
        REAL*8, DIMENSION(:, :),allocatable :: yout
        REAL*8, DIMENSION(:), ALLOCATABLE :: t_eval
        REAL*8, DIMENSION(1) :: params
        INTEGER :: i, midpoint, nparams, nvars
        character(len=100) :: system_name
        params(1) = 0.0d0 ! unused
        ! Initialize arrays
        r = 0.0d0
        W = 0.0d0
        dWdr = 0.0d0
        ! Calculate the break radius
        rbreak = king_break_radius(W0)
        ! First half: solve for W given r
        midpoint = npoints / 2
        ! allocate the size of the output arrays
        nparams = 1
        nvars = 2        
        ALLOCATE(t_eval(midpoint))
        ALLOCATE(yout(nvars,midpoint))
        CALL linspace(0.0d0, rbreak, midpoint, t_eval)
        ! set the initial condition, which is the central potential and zero force at the center
        y0(1) = W0
        y0(2) = 0.0d0
        ! set the size of things
        system_name="king_ode_in_r"
        t_span = [0.0d0, rbreak]
        CALL rk4(system_name, t_span, y0, nparams, params, midpoint, nvars, t_eval, yout)
            ! rk4(system_name,t_span,y0,nparams,params,npoints,nvars,tout,yout)
        r(1:midpoint) = t_eval
        W(1:midpoint) = yout(1,:)
        dWdr(1:midpoint) = yout(2,:)

        ! Second half: solve for r given W
        ! allocate the size of the output arrays
        DEALLOCATE(t_eval)
        DEALLOCATE(yout)

        ALLOCATE(t_eval(npoints - midpoint ))
        ALLOCATE(yout(nvars,npoints - midpoint))
        CALL linspace(W(midpoint),0.0d0, npoints - midpoint , t_eval) ! independent variable is now W

        ! store the rest of the output for W 
        W(midpoint:) = t_eval
        ! set the initial condition, which is the central potential and zero force at the center
        y0(1) = r(midpoint)
        y0(2) = dWdr(midpoint)
        ! set the size of things
        system_name="king_ode_in_w"
        t_span = [W(midpoint), 0.0d0]

        CALL rk4(system_name, t_span, y0, nparams, params, npoints - midpoint , nvars, t_eval, yout)
        ! store the rest of the output for r
        r(1+midpoint:npoints+1) = yout(1,:)
        dWdr(1+midpoint:npoints+1) = yout(2,:)

    END SUBROUTINE solve_king_potential_profile



    FUNCTION KingDensityW(W) RESULT(rho)
        ! King density profile as a function of W
        ! W is the minimium of the potential divided by the velocity dispersion squared
        REAL*8, INTENT(IN) :: W 
        REAL*8 :: rho, sqW, expW, erfW
        REAL*8 :: term1, term2, term3
        sqW = SQRT(W)
        expW = EXP(W)
        erfW = erf_custom(sqW)
        term1 = expW * erfW
        term2 = (2.0/SQRT(PI))*sqW
        term3 = 1.0 + 2.0*W/3.0
        rho = term1 - term2 * term3

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

    SUBROUTINE king_ode_in_r(t,y,dydt,params)
        REAL*8, INTENT(IN) :: t
        REAL*8, INTENT(IN),dimension(2) :: y
        REAL*8, INTENT(OUT),dimension(2) :: dydt
        REAL*8, intent(in), dimension(1) :: params ! unused
        REAL*8 :: r,W, dwdr, d2wdr2, rho 
        ! unpack the variables 
        r = t
        W = y(1)
        dwdr=y(2)
        rho = KingDensityW(W)
        ! numerical gaurd 
        if (r.le.0.0) then 
            d2wdr2=0.0
        else 
            d2wdr2 = -4.0 * PI * rho - 2*dwdr/r
        end if 
        dydt(1) = dwdr
        dydt(2) = d2wdr2

    END SUBROUTINE king_ode_in_r

    SUBROUTINE king_ode_in_w(t,y,dydt,params)
        REAL*8, INTENT(IN) :: t
        REAL*8, INTENT(IN),dimension(2) :: y
        REAL*8, intent(IN), dimension(1) :: params ! unused
        REAL*8, INTENT(OUT),dimension(2) :: dydt
        REAL*8 :: r,W, dwdr, d2wdr2, rho, term1, term2
        ! unpack the variables, in this case W is the independent variable 
        W = t
        r = y(1)
        dwdr=y(2)
        rho = KingDensityW(W)
        term1 = -1.0*(1.0/dwdr)
        term2 = (4.0 * PI * rho + 2.0*dwdr/r)
        d2wdr2 = term1 * term2
        dydt(1) = 1.0/dwdr
        dydt(2) = d2wdr2

    END SUBROUTINE king_ode_in_w


    SUBROUTINE rk4(system_name, t_span, y0, nparams, params, npoints, nvars, tout, yout)
        IMPLICIT NONE
        character*100, intent(in) :: system_name
        REAL*8, INTENT(IN), DIMENSION(2) :: t_span
        REAL*8, INTENT(IN), DIMENSION(nvars) :: y0
        INTEGER, INTENT(IN) :: npoints, nparams, nvars
        REAL*8, INTENT(OUT), DIMENSION(npoints) :: tout
        REAL*8, INTENT(OUT), DIMENSION(nvars, npoints) :: yout
        REAL*8, DIMENSION(nparams), INTENT(IN) :: params
        REAL*8 :: dt, t0, tf
        REAL*8, DIMENSION(nvars) :: k1, k2, k3, k4, y_temp, kterms
        INTEGER :: i
    
        ! Set the initial condition
        t0 = t_span(1)
        tf = t_span(2)
        dt = (tf - t0) / (npoints - 1)
        yout(:, 1) = y0
        tout(1) = t0
    
        if (system_name.eq."exponential_decay") then 
            my_system=>exponential_decay
        else if (system_name.eq."double_system_of_equations") then 
            my_system=>double_system_of_equations
        else if (system_name.eq."king_ode_in_r") then 
            my_system=>king_ode_in_r
        else if (system_name.eq."king_ode_in_w") then 
            my_system=>king_ode_in_w
        else
            print*, system_name, "error, system of equations not implemented"
            stop 
        end if 
        i=1

        ! Runge-Kutta 4th order method
        DO i = 2, npoints
            CALL my_system(tout(i-1), yout(:, i-1), k1, params)
            y_temp = yout(:, i-1) + 0.5d0 * dt * k1
            
            CALL my_system(tout(i-1) + 0.5d0 * dt, y_temp, k2, params)
            y_temp = yout(:, i-1) + 0.5d0 * dt * k2
            
            CALL my_system(tout(i-1) + 0.5d0 * dt, y_temp, k3, params)
            y_temp = yout(:, i-1) + dt * k3
            
            CALL my_system(tout(i-1) + dt, y_temp, k4, params)
            
            kterms = (k1 + (2.0d0 * k2) + (2.0d0 * k3) + k4) / 6.0d0
            yout(:, i) = yout(:, i-1) + kterms * dt
            tout(i) = t0 + (i-1)*dt


        END DO
    END SUBROUTINE rk4

    
    SUBROUTINE linspace(start, end, n, result)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: start, end
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(OUT), DIMENSION(n) :: result
        INTEGER :: i
        REAL*8 :: step
    
        IF (n > 1) THEN
            step = (end - start) / (n - 1)
            DO i = 1, n
                result(i) = start + (i - 1) * step
            END DO
        ELSE
            result(1) = start
        END IF
    END SUBROUTINE linspace
        
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
    