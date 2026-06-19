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

    LOGICAL, PUBLIC :: FLYINGSPHERES_REGISTERED = .FALSE.
    TYPE(flyingsphere_t), ALLOCATABLE, PUBLIC :: flyingsphere_registry(:)
    INTEGER, PUBLIC :: Nflyingspheres = 0 

    CONTAINS

    SUBROUTINE CLEAR()
        FLYINGSPHERES_REGISTERED = .FALSE.
        IF (ALLOCATED(flyingsphere_registry)) DEALLOCATE(flyingsphere_registry)
        Nflyingspheres = 0
    END SUBROUTINE CLEAR

    SUBROUTINE initialize_flyingspheres(n)
        INTEGER, INTENT(IN) :: n
        IF (ALLOCATED(flyingsphere_registry)) DEALLOCATE(flyingsphere_registry)
        ALLOCATE(flyingsphere_registry(n))
        FLYINGSPHERES_REGISTERED = .TRUE.
        Nflyingspheres = n
    END SUBROUTINE initialize_flyingspheres

    SUBROUTINE set_flyingsphere(i,model_name, nparams, structural_parameters, ntimestamps, txyz)
        INTEGER, INTENT(IN) :: i, nparams, ntimestamps
        character(len=64), INTENT(IN):: model_name
        REAl*8, DIMENSION(nparams), INTENT(IN) :: structural_parameters
        REAL*8, DIMENSION(4,ntimestamps),INTENT(IN) :: txyz
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
            print*, "physics spheres no registered"
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

	! REAL*8 FUNCTION force_eval(self,x,y,z)
	! 	CLASS(flyingsphere_t), INTENT(INOUT) :: self
	! 	REAL*8, intent(in) :: x,y,z

	! END FUNCTION force_eval 

end module flyingspheres

! MODULE flyingspheres
! 	! API sketch for a module of orbiting spherical flyingspheres.
! 	!
! 	! Intended shape:
! 	! - many independent flyingsphere objects
! 	! - each object owns its own kinematics table
! 	! - each object may own its own structural parameter tables
! 	! - each object selects a spherical force/potential model
! 	! - the module evaluates all objects at a query time and sums the force
! 	!
! 	! Suggested usage flow:
! 	!   1. clear()
! 	!   2. add_flyingsphere(...)
! 	!   3. configure_flyingsphere_kinematics(...)
! 	!   4. configure_flyingsphere_structure(...)
! 	!   5. optionally configure per-parameter tables
! 	!   6. finalize()
! 	!   7. update_state(t)
! 	!   8. eval_force(...)
! 	!
! 	! Pseudocode for the internal object model:
! 	!
! 	!   type flyingsphere_t
! 	!       character(len=64) :: model_name
! 	!       real*8, allocatable :: t(:)
! 	!       real*8, allocatable :: x(:), y(:), z(:)
! 	!       real*8, allocatable :: vx(:), vy(:), vz(:)
! 	!       real*8, allocatable :: constant_params(:)
! 	!       ! optional per-parameter time tables
! 	!       ! optional cached interpolation indices
! 	!       ! force evaluator callback or model dispatch key
! 	!   end type flyingsphere_t
! 	!
! 	!   type(flyingsphere_t), allocatable :: flyingspheres(:)
! 	!   integer :: Nflyingspheres
! 	!
! 	! Pseudocode for the runtime contract:
! 	!
! 	!   subroutine update_state(t)
! 	!       ! for each flyingsphere:
! 	!       !   interpolate current position from its kinematics table
! 	!       !   interpolate any evolving structural parameters
! 	!       !   store the current state for fast force evaluation
! 	!   end subroutine update_state
! 	!
! 	!   subroutine eval_force(n, x, y, z, ax, ay, az, phi)
! 	!       ! zero accumulators
! 	!       ! for each flyingsphere:
! 	!       !   compute dx, dy, dz relative to the flyingsphere's current position
! 	!       !   dispatch to the selected spherical profile model
! 	!       !   accumulate acceleration and potential
! 	!   end subroutine eval_force
! 	!
! 	! Model examples for the first implementation:
! 	! - plummer
! 	! - hernquist
! 	!
! 	! Future extensions this layout should leave room for:
! 	! - close-encounter diagnostics
! 	! - per-flyingsphere encounter counters
! 	! - more profile types (e.g. NFW, truncated halo, tabulated profiles)
! 	! - time-dependent mass/scale-radius evolution via interpolation tables
	
! 	! USE mathutils, only : linear_interp_scalar, is_strictly_monotonic, is_strictly_decreasing, bracketed_index_search

! 	IMPLICIT NONE 
	
! 	LOGICAL, PUBLIC :: FLYINGSPHERES_REGISTERED = .FALSE.


! 	CONTAINS

! 	SUBROUTINE CLEAR()
! 		FLYINGSPHERES_REGISTERED = .FALSE.

! 	END SUBROUTINE CLEAR

! END MODULE flyingspheres
