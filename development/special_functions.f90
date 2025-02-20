! comment to test git
MODULE special_numerical_functions
    USE, INTRINSIC :: ISO_C_BINDING    
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE f(t, y_in, dydt) BIND(C)
            IMPORT C_DOUBLE
            REAL(C_DOUBLE), INTENT(IN) :: t
            REAL(C_DOUBLE), INTENT(IN), DIMENSION(:) :: y_in
            REAL(C_DOUBLE), INTENT(OUT), DIMENSION(:) :: dydt
        END SUBROUTINE f
    END INTERFACE
    
    CONTAINS

    SUBROUTINE runge_kutta_4_system(y0, t0, tf, h, time, y)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: y0(:), t0, tf, h
        REAL*8, INTENT(OUT), DIMENSION(:), POINTER :: time
        REAL*8, INTENT(OUT), DIMENSION(:,:), POINTER :: y
        INTEGER, INTENT(OUT) :: n
        REAL*8, DIMENSION(:), ALLOCATABLE :: k1, k2, k3, k4, y_temp
        INTEGER :: i, m

        n = INT((tf - t0) / h) + 1
        m = SIZE(y0)

        ALLOCATE(time(n), y(m, n), k1(m), k2(m), k3(m), k4(m), y_temp(m))

        time(1) = t0
        y(:,1) = y0

        DO i = 2, n
            CALL f(time(i-1), y(:, i-1), k1)
            CALL f(time(i-1) + 0.5d0 * h, y(:, i-1) + 0.5d0 * h * k1, k2)
            CALL f(time(i-1) + 0.5d0 * h, y(:, i-1) + 0.5d0 * h * k2, k3)
            CALL f(time(i-1) + h, y(:, i-1) + h * k3, k4)
            y(:, i) = y(:, i-1) + (k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4) * h / 6.0d0
            time(i) = time(i-1) + h
        END DO

        DEALLOCATE(k1, k2, k3, k4, y_temp)
    END SUBROUTINE runge_kutta_4_system

    FUNCTION erf_custom(x) RESULT(y)
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

END MODULE special_numerical_functions

