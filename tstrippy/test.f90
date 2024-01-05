module test
    implicit none
    PUBLIC :: add
    REAL*8, PUBLIC :: result
contains
    SUBROUTINE add(a,b)
        implicit none
        ! local variables
        REAL*8, INTENT(IN) :: a,b
        result = a + b
    END SUBROUTINE add
end module test

