module particles_module
    IMPLICIT NONE
    type :: particles_t 
        integer :: n 
        real*8, allocatable :: x(:,:)
        real*8, allocatable :: v(:,:)
    contains 
        procedure :: allocate
    end type particles_t

    contains 
        subroutine allocate(self, n)
            class(particles_t), intent(inout) :: self
            integer, intent(in) :: n 
            self%n = n 
            allocate(self%x(3,n))
            allocate(self%v(3,n))
        end subroutine allocate
end module particles_module