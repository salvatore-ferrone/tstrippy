MODULE mathutils
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: linear_interp_scalar
    PUBLIC :: bilinear_interp_2d
    PUBLIC :: bicubic_hermite_eval_2d
    PUBLIC :: legendre_p_all_axisymmetric
    PUBLIC :: legendre_dp_dmu_all_axisymmetric
    PUBLIC :: legendre_axisymmetric_basis
    PUBLIC :: gauss_legendre_nodes_weights
    PUBLIC :: bessel_j0_scalar
    PUBLIC :: bessel_j1_scalar

    CONTAINS

    FUNCTION linear_interp_scalar(y0, y1, alpha) RESULT(y)
        REAL*8, INTENT(IN) :: y0, y1, alpha
        REAL*8 :: y
        REAL*8 :: a

        ! Clamp alpha to avoid accidental extrapolation from caller roundoff.
        a = MAX(0.0D0, MIN(1.0D0, alpha))
        y = (1.0D0 - a) * y0 + a * y1
    END FUNCTION linear_interp_scalar

    FUNCTION bilinear_interp_2d(f00, f10, f01, f11, alphaR, alphaZ) RESULT(val)
        ! Bilinear interpolation given four corner values and two fractional
        ! offsets. alphaR and alphaZ are each in [0,1]; the caller is responsible
        ! for clamping them before calling this function.
        !
        !  f01 ------ f11
        !   |          |
        !  f00 ------ f10
        !      alphaR ->
        REAL*8, INTENT(IN) :: f00, f10, f01, f11, alphaR, alphaZ
        REAL*8 :: val
        val = (1.0D0 - alphaR) * (1.0D0 - alphaZ) * f00 &
            +          alphaR  * (1.0D0 - alphaZ) * f10 &
            + (1.0D0 - alphaR) *          alphaZ  * f01 &
            +          alphaR  *          alphaZ  * f11
    END FUNCTION bilinear_interp_2d

    SUBROUTINE bicubic_hermite_eval_2d(f00, f10, f01, f11, &
                                       fx00, fx10, fx01, fx11, &
                                       fy00, fy10, fy01, fy11, &
                                       fxy00, fxy10, fxy01, fxy11, &
                                       dx, dy, tx, ty, f, dfdx, dfdy)
        ! Bicubic Hermite patch evaluation on one rectangle.
        ! Inputs:
        ! - f.. corner function values
        ! - fx.., fy.. corner partial derivatives w.r.t physical x and y
        ! - fxy.. corner mixed derivatives w.r.t physical x then y
        ! - dx, dy rectangle widths in physical coordinates
        ! - tx, ty local coordinates in [0,1]
        ! Outputs:
        ! - f, dfdx, dfdy evaluated from one single interpolant.
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: f00, f10, f01, f11
        REAL*8, INTENT(IN) :: fx00, fx10, fx01, fx11
        REAL*8, INTENT(IN) :: fy00, fy10, fy01, fy11
        REAL*8, INTENT(IN) :: fxy00, fxy10, fxy01, fxy11
        REAL*8, INTENT(IN) :: dx, dy, tx, ty
        REAL*8, INTENT(OUT) :: f, dfdx, dfdy

        REAL*8 :: t, u
        REAL*8 :: h00_t, h10_t, h01_t, h11_t
        REAL*8 :: h00_u, h10_u, h01_u, h11_u
        REAL*8 :: dh00_t, dh10_t, dh01_t, dh11_t
        REAL*8 :: dh00_u, dh10_u, dh01_u, dh11_u
        REAL*8 :: f0, f1, fy0, fy1
        REAL*8 :: dfdx0, dfdx1, dfdxy0, dfdxy1

        t = MAX(0.0D0, MIN(1.0D0, tx))
        u = MAX(0.0D0, MIN(1.0D0, ty))

        ! Cubic Hermite basis in x (t)
        h00_t =  2.0D0*t**3 - 3.0D0*t**2 + 1.0D0
        h10_t =        t**3 - 2.0D0*t**2 + t
        h01_t = -2.0D0*t**3 + 3.0D0*t**2
        h11_t =        t**3 -       t**2

        dh00_t =  6.0D0*t**2 - 6.0D0*t
        dh10_t =  3.0D0*t**2 - 4.0D0*t + 1.0D0
        dh01_t = -6.0D0*t**2 + 6.0D0*t
        dh11_t =  3.0D0*t**2 - 2.0D0*t

        ! Cubic Hermite basis in y (u)
        h00_u =  2.0D0*u**3 - 3.0D0*u**2 + 1.0D0
        h10_u =        u**3 - 2.0D0*u**2 + u
        h01_u = -2.0D0*u**3 + 3.0D0*u**2
        h11_u =        u**3 -       u**2

        dh00_u =  6.0D0*u**2 - 6.0D0*u
        dh10_u =  3.0D0*u**2 - 4.0D0*u + 1.0D0
        dh01_u = -6.0D0*u**2 + 6.0D0*u
        dh11_u =  3.0D0*u**2 - 2.0D0*u

        ! Interpolate in x at y=0 and y=1, including y-derivative tracks.
        f0  = h00_t*f00  + h10_t*dx*fx00  + h01_t*f10  + h11_t*dx*fx10
        f1  = h00_t*f01  + h10_t*dx*fx01  + h01_t*f11  + h11_t*dx*fx11
        fy0 = h00_t*fy00 + h10_t*dx*fxy00 + h01_t*fy10 + h11_t*dx*fxy10
        fy1 = h00_t*fy01 + h10_t*dx*fxy01 + h01_t*fy11 + h11_t*dx*fxy11

        ! d/dx of the x-interpolants.
        dfdx0  = (dh00_t*f00  + dh10_t*dx*fx00  + dh01_t*f10  + dh11_t*dx*fx10) / dx
        dfdx1  = (dh00_t*f01  + dh10_t*dx*fx01  + dh01_t*f11  + dh11_t*dx*fx11) / dx
        dfdxy0 = (dh00_t*fy00 + dh10_t*dx*fxy00 + dh01_t*fy10 + dh11_t*dx*fxy10) / dx
        dfdxy1 = (dh00_t*fy01 + dh10_t*dx*fxy01 + dh01_t*fy11 + dh11_t*dx*fxy11) / dx

        ! Final Hermite interpolation in y.
        f    = h00_u*f0    + h10_u*dy*fy0    + h01_u*f1    + h11_u*dy*fy1
        dfdx = h00_u*dfdx0 + h10_u*dy*dfdxy0 + h01_u*dfdx1 + h11_u*dy*dfdxy1
        dfdy = (dh00_u*f0   + dh10_u*dy*fy0   + dh01_u*f1   + dh11_u*dy*fy1) / dy
    END SUBROUTINE bicubic_hermite_eval_2d

    SUBROUTINE legendre_p_all_axisymmetric(lmax, mu, p)
        ! Compute P_l(mu) for l=0..lmax using three-term recurrence.
        INTEGER, INTENT(IN) :: lmax
        REAL*8, INTENT(IN) :: mu
        REAL*8, INTENT(OUT), DIMENSION(0:lmax) :: p
        INTEGER :: l

        p(0) = 1.0D0
        IF (lmax == 0) RETURN

        p(1) = mu
        IF (lmax == 1) RETURN

        DO l = 1, lmax - 1
            p(l+1) = ((2.0D0*l + 1.0D0) * mu * p(l) - l * p(l-1)) / (l + 1.0D0)
        END DO
    END SUBROUTINE legendre_p_all_axisymmetric

    SUBROUTINE legendre_dp_dmu_all_axisymmetric(lmax, mu, p, dp)
        ! Compute dP_l/dmu for l=0..lmax from P_l values.
        ! Uses analytic endpoint limits at mu=+/-1 to avoid singular forms.
        INTEGER, INTENT(IN) :: lmax
        REAL*8, INTENT(IN) :: mu
        REAL*8, INTENT(IN), DIMENSION(0:lmax) :: p
        REAL*8, INTENT(OUT), DIMENSION(0:lmax) :: dp
        INTEGER :: l
        REAL*8 :: den
        REAL*8, PARAMETER :: eps_mu = 1.0D-12

        dp(0) = 0.0D0
        IF (lmax == 0) RETURN

        IF (ABS(mu - 1.0D0) < eps_mu) THEN
            DO l = 1, lmax
                dp(l) = 0.5D0 * l * (l + 1.0D0)
            END DO
            RETURN
        END IF

        IF (ABS(mu + 1.0D0) < eps_mu) THEN
            DO l = 1, lmax
                IF (MOD(l + 1, 2) == 0) THEN
                    dp(l) = 0.5D0 * l * (l + 1.0D0)
                ELSE
                    dp(l) = -0.5D0 * l * (l + 1.0D0)
                END IF
            END DO
            RETURN
        END IF

        den = mu*mu - 1.0D0
        DO l = 1, lmax
            dp(l) = l * (mu * p(l) - p(l-1)) / den
        END DO
    END SUBROUTINE legendre_dp_dmu_all_axisymmetric

    SUBROUTINE legendre_axisymmetric_basis(lmax, mu, p, dp)
        ! Convenience wrapper to compute both P_l and dP_l/dmu for l=0..lmax.
        INTEGER, INTENT(IN) :: lmax
        REAL*8, INTENT(IN) :: mu
        REAL*8, INTENT(OUT), DIMENSION(0:lmax) :: p, dp

        CALL legendre_p_all_axisymmetric(lmax, mu, p)
        CALL legendre_dp_dmu_all_axisymmetric(lmax, mu, p, dp)
    END SUBROUTINE legendre_axisymmetric_basis

    SUBROUTINE gauss_legendre_nodes_weights(n, nodes, weights)
        ! Gauss-Legendre nodes and weights on [-1, 1] via Newton-Raphson
        ! on roots of P_n. Nodes are returned in ascending order.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(OUT), DIMENSION(n) :: nodes, weights
        REAL*8, PARAMETER :: pi_gl = 3.14159265358979323846D0
        REAL*8, PARAMETER :: eps_gl = 1.0D-15
        INTEGER :: i, j, m
        REAL*8 :: x, x_old, p_prev, p_curr, p_next, dp

        m = (n + 1) / 2
        DO i = 1, m
            x = COS(pi_gl * (i - 0.25D0) / (n + 0.5D0))
            DO
                p_prev = 1.0D0
                p_curr = x
                DO j = 2, n
                    p_next = ((2*j - 1) * x * p_curr - (j - 1) * p_prev) / DBLE(j)
                    p_prev = p_curr
                    p_curr = p_next
                END DO
                dp = n * (x * p_curr - p_prev) / (x*x - 1.0D0)
                x_old = x
                x = x_old - p_curr / dp
                IF (ABS(x - x_old) <= eps_gl) EXIT
            END DO
            nodes(i)       = -x
            nodes(n+1-i)   =  x
            weights(i)     = 2.0D0 / ((1.0D0 - x*x) * dp*dp)
            weights(n+1-i) = weights(i)
        END DO
    END SUBROUTINE gauss_legendre_nodes_weights

    FUNCTION bessel_j0_scalar(x) RESULT(y)
        ! Wrapper around compiler intrinsic BESSEL_J0.
        REAL*8, INTENT(IN) :: x
        REAL*8 :: y
        y = BESSEL_J0(x)
    END FUNCTION bessel_j0_scalar

    FUNCTION bessel_j1_scalar(x) RESULT(y)
        ! Wrapper around compiler intrinsic BESSEL_J1.
        REAL*8, INTENT(IN) :: x
        REAL*8 :: y
        y = BESSEL_J1(x)
    END FUNCTION bessel_j1_scalar

END MODULE mathutils
