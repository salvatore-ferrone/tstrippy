! the god module for orchestating simulations

module simulator_module
    use particles_module
    IMPLICIT NONE
    private 
    public :: create
    integer, parameter ::  MAX_SIMS = 8 
    type :: simulator_t
        class(particles_t), allocatable :: particles 
    CONTAINS
        procedure :: set_initial_conditions
    end type simulator_t 


    type(simulator_t), save :: SIM 
    logical, save :: SIM_INITIALIZE = .FALSE.

    CONTAINS

        ! PUBLIC FACING API 
        subroutine create()
            print*, "create"
        end subroutine create 

        ! INTENRAL 
        subroutine set_initial_conditions(n,x,v)
            INTEGER, INTENT(IN) :: n 
            REAL*8, INTENT(IN), DIMENSION(3,n) :: x,v
            print*, "hello"
            SIM%particles%allocat
        end subroutine set_initial_conditions


end module simulator_module