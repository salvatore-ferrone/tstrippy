MODULE mathutils
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: linear_interp_scalar

    CONTAINS

    FUNCTION linear_interp_scalar(y0, y1, alpha) RESULT(y)
        REAL*8, INTENT(IN) :: y0, y1, alpha
        REAL*8 :: y
        REAL*8 :: a

        ! Clamp alpha to avoid accidental extrapolation from caller roundoff.
        a = MAX(0.0D0, MIN(1.0D0, alpha))
        y = (1.0D0 - a) * y0 + a * y1
    END FUNCTION linear_interp_scalar

END MODULE mathutils
