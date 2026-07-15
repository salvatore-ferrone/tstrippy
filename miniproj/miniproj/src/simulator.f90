! the god module for orchestating simulations

module simulator_module
    use particles_module
    IMPLICIT NONE

    type :: simulator_t
        class(particles_t), allocatable :: particles 
    CONTAINS
        procedure :: setinitialconditions
    end type simulator_t 

    CONTAINS

        subroutine setinitialconditions(self)
            class(simulator_t), intent(in) :: self 
            print*, "hello"
        end subroutine setinitialconditions

end module simulator_module