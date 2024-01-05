module integrator
    implicit none
    public :: addition
contains
    SUBROUTINE addition(a,b,c)
        real, intent(in) :: a,b
        real, intent(out) :: c
        c = a + b
    END SUBROUTINE addition
end module integrator