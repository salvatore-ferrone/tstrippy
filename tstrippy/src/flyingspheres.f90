MODULE flyingspheres
	! API sketch for a module of orbiting spherical perturbers.
	!
	! Intended shape:
	! - many independent perturber objects
	! - each object owns its own kinematics table
	! - each object may own its own structural parameter tables
	! - each object selects a spherical force/potential model
	! - the module evaluates all objects at a query time and sums the force
	!
	! Suggested usage flow:
	!   1. clear()
	!   2. add_perturber(...)
	!   3. configure_perturber_kinematics(...)
	!   4. configure_perturber_structure(...)
	!   5. optionally configure per-parameter tables
	!   6. finalize()
	!   7. update_state(t)
	!   8. eval_force(...)
	!
	! Pseudocode for the internal object model:
	!
	!   type perturber_t
	!       character(len=64) :: model_name
	!       real*8, allocatable :: t(:)
	!       real*8, allocatable :: x(:), y(:), z(:)
	!       real*8, allocatable :: vx(:), vy(:), vz(:)
	!       real*8, allocatable :: constant_params(:)
	!       ! optional per-parameter time tables
	!       ! optional cached interpolation indices
	!       ! force evaluator callback or model dispatch key
	!   end type perturber_t
	!
	!   type(perturber_t), allocatable :: perturbers(:)
	!   integer :: nperturbers
	!
	! Pseudocode for the runtime contract:
	!
	!   subroutine update_state(t)
	!       ! for each perturber:
	!       !   interpolate current position from its kinematics table
	!       !   interpolate any evolving structural parameters
	!       !   store the current state for fast force evaluation
	!   end subroutine update_state
	!
	!   subroutine eval_force(n, x, y, z, ax, ay, az, phi)
	!       ! zero accumulators
	!       ! for each perturber:
	!       !   compute dx, dy, dz relative to the perturber's current position
	!       !   dispatch to the selected spherical profile model
	!       !   accumulate acceleration and potential
	!   end subroutine eval_force
	!
	! Model examples for the first implementation:
	! - plummer
	! - hernquist
	!
	! Future extensions this layout should leave room for:
	! - close-encounter diagnostics
	! - per-perturber encounter counters
	! - more profile types (e.g. NFW, truncated halo, tabulated profiles)
	! - time-dependent mass/scale-radius evolution via interpolation tables

END MODULE flyingspheres
