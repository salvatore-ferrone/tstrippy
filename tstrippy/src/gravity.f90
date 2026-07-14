module gravity
    implicit none
    
contains
    subroutine add_component(model_name, params, backend, config)
        !f2py intent(in) :: model_name
        !f2py intent(in) :: params
        !f2py optional :: backend
        !f2py optional :: config

        character(len=*), intent(in) :: model_name
        real(8), intent(in) :: params(:)

        character(len=*), optional, intent(in) :: backend
        character(len=*), optional, intent(in) :: config

        integer :: i

        print*, "model_name", model_name

        do i=1,size(params)
            print*, "params(i)", i, params(i)
        end do

        if (present(backend)) then
            print *, "backend is present"
            print*, backend
        end if 

        if (present(config)) then
            print *, "config is present"
            print*, config
        end if 
    end subroutine add_component


end module gravity