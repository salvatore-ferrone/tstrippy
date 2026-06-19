MODULE flyingspheres
    USE mathutils, only : linear_interp_scalar, bracketed_index_search
    IMPLICIT NONE 

	! the function headers
    ABSTRACT INTERFACE
        SUBROUTINE force_eval_iface(params,n,x,y,z,force)
            REAL*8, INTENT(IN) :: params(:)
            INTEGER, INTENT(IN) :: n
            REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
            REAL*8, INTENT(OUT), DIMENSION(n,3) :: force
        END SUBROUTINE force_eval_iface
    END INTERFACE

    type, private :: flyingsphere_t
        CHARACTER(len=64) :: model_name
        REAL*8, ALLOCATABLE :: txyz(:,:)
        REAL*8, ALLOCATABLE :: structural_parameters(:)
        REAL*8 :: t,x,y,z ! the current values
        integer :: index = 1 
        PROCEDURE(force_eval_iface), pointer, NOPASS, private :: model_force => NULL()
    end type flyingsphere_t

    type, private :: parameter_table_t 
        REAL*8 :: value = 0.0D0
        INTEGER :: index = 1 
        REAL*8, ALLOCATABLE :: t(:), values(:)
        contains 
            PROCEDURE :: update_parameter
    end type parameter_table_t


    REAL*8, PARAMETER, PRIVATE :: G_DEFAULT = 4.30091727D-6
    REAL*8, PUBLIC  :: G_FLYINGSPHERES = G_DEFAULT
    LOGICAL, PUBLIC :: G_IS_DEFAULT = .TRUE.

    LOGICAL, PUBLIC :: FLYINGSPHERES_REGISTERED = .FALSE.
    LOGICAL, PUBLIC :: FLYINGSPHERES_FINALIZED = .FALSE.
    TYPE(flyingsphere_t), ALLOCATABLE, PUBLIC :: flyingsphere_registry(:)
    INTEGER, PUBLIC :: Nflyingspheres = 0 

    CONTAINS

    SUBROUTINE CLEAR()
        FLYINGSPHERES_REGISTERED = .FALSE.
        FLYINGSPHERES_FINALIZED = .FALSE.
        IF (ALLOCATED(flyingsphere_registry)) DEALLOCATE(flyingsphere_registry)
        Nflyingspheres = 0
    END SUBROUTINE CLEAR

    SUBROUTINE finalize_flyingspheres()
        INTEGER :: i
        FLYINGSPHERES_FINALIZED = .FALSE.

        IF (.NOT. FLYINGSPHERES_REGISTERED) THEN
            PRINT*, "ERROR: finalize_flyingspheres called before initialize_flyingspheres"
            RETURN
        END IF

        IF (.NOT. ALLOCATED(flyingsphere_registry)) THEN
            PRINT*, "ERROR: finalize_flyingspheres: registry is not allocated"
            RETURN
        END IF

        IF (Nflyingspheres < 1) THEN
            PRINT*, "ERROR: finalize_flyingspheres: no flyingspheres configured"
            RETURN
        END IF

        DO i = 1, Nflyingspheres
            IF (.NOT. ALLOCATED(flyingsphere_registry(i)%structural_parameters)) THEN
                PRINT*, "ERROR: finalize_flyingspheres: structural_parameters missing for index", i
                RETURN
            END IF
            IF (.NOT. ALLOCATED(flyingsphere_registry(i)%txyz)) THEN
                PRINT*, "ERROR: finalize_flyingspheres: txyz missing for index", i
                RETURN
            END IF
            IF (SIZE(flyingsphere_registry(i)%txyz, 1) /= 4) THEN
                PRINT*, "ERROR: finalize_flyingspheres: txyz must have shape (4,ntimestamps) for index", i
                RETURN
            END IF
            IF (SIZE(flyingsphere_registry(i)%txyz, 2) < 2) THEN
                PRINT*, "ERROR: finalize_flyingspheres: ntimestamps must be >= 2 for index", i
                RETURN
            END IF
            SELECT CASE (TRIM(flyingsphere_registry(i)%model_name))
            CASE ("plummer")
                flyingsphere_registry(i)%model_force => plummer_force
            ! CASE ("hernquist")
                ! flyingsphere_registry(i)%model_force => hernquist_force
            CASE DEFAULT
                PRINT*, "ERROR: unknown flyingsphere model: ", TRIM(flyingsphere_registry(i)%model_name)
                RETURN
            END SELECT
        END DO
        FLYINGSPHERES_FINALIZED = .TRUE.
    END SUBROUTINE finalize_flyingspheres

    SUBROUTINE initialize_flyingspheres(n)
        INTEGER, INTENT(IN) :: n
        IF (ALLOCATED(flyingsphere_registry)) DEALLOCATE(flyingsphere_registry)
        ALLOCATE(flyingsphere_registry(n))
        FLYINGSPHERES_REGISTERED = .TRUE.
        FLYINGSPHERES_FINALIZED = .FALSE.
        Nflyingspheres = n
    END SUBROUTINE initialize_flyingspheres

    SUBROUTINE set_flyingsphere(i,model_name, nparams, structural_parameters, ntimestamps, txyz)
        INTEGER, INTENT(IN) :: i, nparams, ntimestamps
        character(len=64), INTENT(IN):: model_name
        REAl*8, DIMENSION(nparams), INTENT(IN) :: structural_parameters
        REAL*8, DIMENSION(4,ntimestamps),INTENT(IN) :: txyz
        FLYINGSPHERES_FINALIZED = .FALSE.
        IF (ALLOCATED(flyingsphere_registry(i)%structural_parameters)) DEALLOCATE(flyingsphere_registry(i)%structural_parameters)
        IF (ALLOCATED(flyingsphere_registry(i)%txyz)) DEALLOCATE(flyingsphere_registry(i)%txyz)
        allocate(flyingsphere_registry(i)%txyz(4,ntimestamps))
        allocate(flyingsphere_registry(i)%structural_parameters(nparams))
        flyingsphere_registry(i)%model_name = TRIM(model_name)
        flyingsphere_registry(i)%structural_parameters = structural_parameters
        flyingsphere_registry(i)%txyz = txyz
    end subroutine set_flyingsphere

    SUBROUTINE update_state(t)
        REAL*8, intent(in) :: t 
        INTEGER :: i
        REAL*8 :: alpha, T0, TF, DT
        
        if (.NOT. FLYINGSPHERES_REGISTERED) THEN 
            print*, "FLYINGSPHERES_REGISTERED is false! cannot update_state"
            return 
        end if 

        DO i=1,Nflyingspheres
            flyingsphere_registry(i)%index=bracketed_index_search(t,flyingsphere_registry(i)%index,flyingsphere_registry(i)%txyz(1,:))
            T0=flyingsphere_registry(i)%txyz(1,flyingsphere_registry(i)%index)
            TF=flyingsphere_registry(i)%txyz(1,flyingsphere_registry(i)%index+1)
            DT = TF-T0
            alpha=t-t0
            flyingsphere_registry(i)%t=t
            flyingsphere_registry(i)%x=linear_interp_scalar(&
                flyingsphere_registry(i)%txyz(2,flyingsphere_registry(i)%index),&
                flyingsphere_registry(i)%txyz(2,flyingsphere_registry(i)%index+1),&
                alpha)
            flyingsphere_registry(i)%y=linear_interp_scalar(&
                flyingsphere_registry(i)%txyz(3,flyingsphere_registry(i)%index),&
                flyingsphere_registry(i)%txyz(3,flyingsphere_registry(i)%index+1),&
                alpha)
            flyingsphere_registry(i)%z=linear_interp_scalar(&
                flyingsphere_registry(i)%txyz(4,flyingsphere_registry(i)%index),&
                flyingsphere_registry(i)%txyz(4,flyingsphere_registry(i)%index+1),&
                alpha)
        end do 
		
		! eventually implement updating specific parameters

    END subroutine update_state

    subroutine update_parameter(self,t)
        CLASS(parameter_table_t), INTENT(INOUT) :: self
        REAL*8 :: t 
        REAL*8 :: alpha, T0, TF, DT
        self%index = bracketed_index_search(t,self%index, self%t)
        T0 = self%t(self%index)
        TF = self%t(self%index+1)
        DT = TF-T0
        alpha=t-t0
        self%value = linear_interp_scalar(self%values(self%index),&
        self%values(self%index+1), alpha)
    END subroutine update_parameter

    SUBROUTINE eval_force_from_one_flyingsphere(i, nparticles, x_particles, y_particles, z_particles, ax, ay, az)
        INTEGER, INTENT(IN) :: i, nparticles
        REAL*8, INTENT(IN) :: x_particles(nparticles), y_particles(nparticles), z_particles(nparticles)
        REAL*8, INTENT(OUT) :: ax(nparticles), ay(nparticles), az(nparticles)
        
        REAL*8 :: dx(nparticles), dy(nparticles), dz(nparticles)
        REAL*8 :: force_temp(nparticles, 3)
        
        ! Relative positions
        dx = x_particles - flyingsphere_registry(i)%x
        dy = y_particles - flyingsphere_registry(i)%y
        dz = z_particles - flyingsphere_registry(i)%z
        
        ! Call this flyingsphere's model with its params and relative positions
        CALL flyingsphere_registry(i)%model_force( &
            flyingsphere_registry(i)%structural_parameters, &
            nparticles, dx, dy, dz, force_temp)
        
        ax = force_temp(:, 1)
        ay = force_temp(:, 2)
        az = force_temp(:, 3)
    END SUBROUTINE

    SUBROUTINE eval_total_force(nparticles, x, y, z, ax, ay, az)
        INTEGER, INTENT(IN) :: nparticles
        REAL*8, INTENT(IN) :: x(nparticles), y(nparticles), z(nparticles)
        REAL*8, INTENT(OUT) :: ax(nparticles), ay(nparticles), az(nparticles)
        
        INTEGER :: i
        REAL*8 :: ax_tmp(nparticles), ay_tmp(nparticles), az_tmp(nparticles)
        
        ax = 0.0D0
        ay = 0.0D0
        az = 0.0D0
        
        DO i = 1, Nflyingspheres
            IF (.NOT. ASSOCIATED(flyingsphere_registry(i)%model_force)) CYCLE
            CALL eval_force_from_one_flyingsphere(i, nparticles, x, y, z, ax_tmp, ay_tmp, az_tmp)
            ax = ax + ax_tmp
            ay = ay + ay_tmp
            az = az + az_tmp
        END DO
    END SUBROUTINE

	! ANALYTICAL MODELS
    SUBROUTINE plummer_force(params, n, x, y, z, force)
        REAL*8, DIMENSION(:), INTENT(IN) :: params
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n,3) :: force
        REAL*8, DIMENSION(n) :: r, amod
        REAL*8 :: m, b

        IF (SIZE(params) < 2) THEN
            force = 0.0D0
            RETURN
        END IF

        m = params(1)
        b = params(2)
        r = SQRT(x*x + y*y + z*z)
        amod = -G_FLYINGSPHERES * m / (r*r + b*b)**1.5

        force(:,1) = amod*x
        force(:,2) = amod*y
        force(:,3) = amod*z
    END SUBROUTINE plummer_force

end module flyingspheres

