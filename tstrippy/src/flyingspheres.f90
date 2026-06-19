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

    type, private :: perturber_t
        CHARACTER(len=64) :: model_name
        REAL*8, ALLOCATABLE :: txyz(:,:)
        REAL*8, ALLOCATABLE :: current_model_params(:)
        PROCEDURE(force_eval_iface), pointer, NOPASS, private :: model_force => NULL()
        ! contains 
            ! PROCEDURE :: force_eval
    end type perturber_t

    type, private :: parameter_table_t 
        REAL*8 :: value = 0.0D0
        INTEGER :: index = 1 
        REAL*8, ALLOCATABLE :: t(:), values(:)
        contains 
            PROCEDURE :: update_parameter
    end type parameter_table_t


    LOGICAL, PUBLIC :: FLYINGSPHERES_REGISTERED = .FALSE.
    TYPE(perturber_t), ALLOCATABLE, PUBLIC :: perturbers(:)
    INTEGER, PUBLIC :: NPERTURBERS = 0 


    CONTAINS

    SUBROUTINE CLEAR()
        FLYINGSPHERES_REGISTERED = .FALSE.
        IF (ALLOCATED(perturbers)) DEALLOCATE(perturbers)
        NPERTURBERS = 0
        print*, "okay"
    END SUBROUTINE CLEAR

    SUBROUTINE INIT_PERTURBERS(n)
        INTEGER, INTENT(IN) :: n
        IF (ALLOCATED(perturbers)) DEALLOCATE(perturbers)
        ALLOCATE(perturbers(n))
        FLYINGSPHERES_REGISTERED = .TRUE.
        NPERTURBERS = n
    END SUBROUTINE INIT_PERTURBERS

	SUBROUTINE set_perturber(i,model_name,structural_parameters, txyz)
		INTEGER, INTENT(IN) :: i 
		character(len=64), INTENT(IN):: model_name
		REAl*8, DIMENSION(:) :: structural_parameters
		REAL*8, DIMENSION(:,:) :: txyz

		print*, "ok"

	end subroutine set_perturber

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
	! 	CLASS(perturber_t), INTENT(INOUT) :: self
	! 	REAL*8, intent(in) :: x,y,z

	! END FUNCTION force_eval 

end module flyingspheres

! MODULE flyingspheres
! 	! API sketch for a module of orbiting spherical perturbers.
! 	!
! 	! Intended shape:
! 	! - many independent perturber objects
! 	! - each object owns its own kinematics table
! 	! - each object may own its own structural parameter tables
! 	! - each object selects a spherical force/potential model
! 	! - the module evaluates all objects at a query time and sums the force
! 	!
! 	! Suggested usage flow:
! 	!   1. clear()
! 	!   2. add_perturber(...)
! 	!   3. configure_perturber_kinematics(...)
! 	!   4. configure_perturber_structure(...)
! 	!   5. optionally configure per-parameter tables
! 	!   6. finalize()
! 	!   7. update_state(t)
! 	!   8. eval_force(...)
! 	!
! 	! Pseudocode for the internal object model:
! 	!
! 	!   type perturber_t
! 	!       character(len=64) :: model_name
! 	!       real*8, allocatable :: t(:)
! 	!       real*8, allocatable :: x(:), y(:), z(:)
! 	!       real*8, allocatable :: vx(:), vy(:), vz(:)
! 	!       real*8, allocatable :: constant_params(:)
! 	!       ! optional per-parameter time tables
! 	!       ! optional cached interpolation indices
! 	!       ! force evaluator callback or model dispatch key
! 	!   end type perturber_t
! 	!
! 	!   type(perturber_t), allocatable :: perturbers(:)
! 	!   integer :: nperturbers
! 	!
! 	! Pseudocode for the runtime contract:
! 	!
! 	!   subroutine update_state(t)
! 	!       ! for each perturber:
! 	!       !   interpolate current position from its kinematics table
! 	!       !   interpolate any evolving structural parameters
! 	!       !   store the current state for fast force evaluation
! 	!   end subroutine update_state
! 	!
! 	!   subroutine eval_force(n, x, y, z, ax, ay, az, phi)
! 	!       ! zero accumulators
! 	!       ! for each perturber:
! 	!       !   compute dx, dy, dz relative to the perturber's current position
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
! 	! - per-perturber encounter counters
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
