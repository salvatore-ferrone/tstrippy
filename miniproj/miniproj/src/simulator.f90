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
        subroutine setinitialconditions(self, x,y,z,vx,vy,vz)
            class(simulator_t), intent(in)
            print*, hello