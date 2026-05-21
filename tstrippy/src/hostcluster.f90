MODULE hostcluster
    IMPLICIT NONE

    INTEGER, PARAMETER, PUBLIC :: HOST_PARAM_MODE_NONE = 0
    INTEGER, PARAMETER, PUBLIC :: HOST_PARAM_MODE_CONSTANT = 1
    INTEGER, PARAMETER, PUBLIC :: HOST_PARAM_MODE_LAW = 2
    INTEGER, PARAMETER, PUBLIC :: HOST_PARAM_MODE_TABLE = 3

    LOGICAL, PUBLIC :: HOST_REGISTERED = .FALSE.
    LOGICAL, PUBLIC :: HOST_FINALIZED = .FALSE.
    LOGICAL, PUBLIC :: HOST_KINEMATICS_SET = .FALSE.
    INTEGER, PUBLIC :: HOST_PARAM_MODE = HOST_PARAM_MODE_NONE

    CHARACTER(LEN=64), PUBLIC :: HOST_BACKEND_NAME = ""
    CHARACTER(LEN=64), PUBLIC :: HOST_LAW_NAME = ""

    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_times
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_x, host_y, host_z
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_vx, host_vy, host_vz

    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_params_constant
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_law_params
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_param_times
    REAL*8, DIMENSION(:,:), PUBLIC, ALLOCATABLE :: host_param_table

    REAL*8, PUBLIC :: host_x_current = 0.0D0
    REAL*8, PUBLIC :: host_y_current = 0.0D0
    REAL*8, PUBLIC :: host_z_current = 0.0D0
    REAL*8, PUBLIC :: host_vx_current = 0.0D0
    REAL*8, PUBLIC :: host_vy_current = 0.0D0
    REAL*8, PUBLIC :: host_vz_current = 0.0D0

    PUBLIC :: clear
    PUBLIC :: add_hostcluster
    PUBLIC :: init_hostcluster_kinematics
    PUBLIC :: set_hostcluster_backend
    PUBLIC :: set_hostcluster_structure_constant
    PUBLIC :: set_hostcluster_structure_law
    PUBLIC :: set_hostcluster_structure_table
    PUBLIC :: finalize_hostcluster
    PUBLIC :: update_hostcluster_state
    PUBLIC :: force_hostcluster_on_particles

CONTAINS

    SUBROUTINE clear()
        IF (ALLOCATED(host_times)) DEALLOCATE(host_times)
        IF (ALLOCATED(host_x)) DEALLOCATE(host_x)
        IF (ALLOCATED(host_y)) DEALLOCATE(host_y)
        IF (ALLOCATED(host_z)) DEALLOCATE(host_z)
        IF (ALLOCATED(host_vx)) DEALLOCATE(host_vx)
        IF (ALLOCATED(host_vy)) DEALLOCATE(host_vy)
        IF (ALLOCATED(host_vz)) DEALLOCATE(host_vz)

        IF (ALLOCATED(host_params_constant)) DEALLOCATE(host_params_constant)
        IF (ALLOCATED(host_law_params)) DEALLOCATE(host_law_params)
        IF (ALLOCATED(host_param_times)) DEALLOCATE(host_param_times)
        IF (ALLOCATED(host_param_table)) DEALLOCATE(host_param_table)

        HOST_REGISTERED = .FALSE.
        HOST_FINALIZED = .FALSE.
        HOST_KINEMATICS_SET = .FALSE.
        HOST_PARAM_MODE = HOST_PARAM_MODE_NONE
        HOST_BACKEND_NAME = ""
        HOST_LAW_NAME = ""
        host_x_current = 0.0D0
        host_y_current = 0.0D0
        host_z_current = 0.0D0
        host_vx_current = 0.0D0
        host_vy_current = 0.0D0
        host_vz_current = 0.0D0
    END SUBROUTINE clear

    SUBROUTINE add_hostcluster()
        HOST_REGISTERED = .TRUE.
        HOST_FINALIZED = .FALSE.
    END SUBROUTINE add_hostcluster

    SUBROUTINE init_hostcluster_kinematics(ntimes, t, x, y, z, vx, vy, vz)
        INTEGER, INTENT(IN) :: ntimes
        REAL*8, INTENT(IN), DIMENSION(ntimes) :: t, x, y, z, vx, vy, vz

        IF (.NOT. HOST_REGISTERED) THEN
            PRINT*, "WARNING: hostcluster not registered. Call add_hostcluster first"
            RETURN
        END IF

        IF (ntimes < 2) THEN
            PRINT*, "WARNING: init_hostcluster_kinematics requires ntimes >= 2"
            RETURN
        END IF

        IF (ALLOCATED(host_times)) DEALLOCATE(host_times)
        IF (ALLOCATED(host_x)) DEALLOCATE(host_x)
        IF (ALLOCATED(host_y)) DEALLOCATE(host_y)
        IF (ALLOCATED(host_z)) DEALLOCATE(host_z)
        IF (ALLOCATED(host_vx)) DEALLOCATE(host_vx)
        IF (ALLOCATED(host_vy)) DEALLOCATE(host_vy)
        IF (ALLOCATED(host_vz)) DEALLOCATE(host_vz)

        ALLOCATE(host_times(ntimes), host_x(ntimes), host_y(ntimes), host_z(ntimes))
        ALLOCATE(host_vx(ntimes), host_vy(ntimes), host_vz(ntimes))

        host_times = t
        host_x = x
        host_y = y
        host_z = z
        host_vx = vx
        host_vy = vy
        host_vz = vz

        host_x_current = x(1)
        host_y_current = y(1)
        host_z_current = z(1)
        host_vx_current = vx(1)
        host_vy_current = vy(1)
        host_vz_current = vz(1)

        HOST_KINEMATICS_SET = .TRUE.
        HOST_FINALIZED = .FALSE.
    END SUBROUTINE init_hostcluster_kinematics

    SUBROUTINE set_hostcluster_backend(model_name)
        CHARACTER(LEN=*), INTENT(IN) :: model_name

        IF (.NOT. HOST_REGISTERED) THEN
            PRINT*, "WARNING: hostcluster not registered. Call add_hostcluster first"
            RETURN
        END IF

        HOST_BACKEND_NAME = TRIM(model_name)
        HOST_FINALIZED = .FALSE.
    END SUBROUTINE set_hostcluster_backend

    SUBROUTINE set_hostcluster_structure_constant(params, nparams)
        INTEGER, INTENT(IN) :: nparams
        REAL*8, INTENT(IN), DIMENSION(nparams) :: params

        IF (.NOT. HOST_REGISTERED) THEN
            PRINT*, "WARNING: hostcluster not registered. Call add_hostcluster first"
            RETURN
        END IF

        IF (nparams < 1) THEN
            PRINT*, "WARNING: set_hostcluster_structure_constant requires nparams >= 1"
            RETURN
        END IF

        IF (ALLOCATED(host_params_constant)) DEALLOCATE(host_params_constant)
        ALLOCATE(host_params_constant(nparams))
        host_params_constant = params

        HOST_PARAM_MODE = HOST_PARAM_MODE_CONSTANT
        HOST_FINALIZED = .FALSE.
    END SUBROUTINE set_hostcluster_structure_constant

    SUBROUTINE set_hostcluster_structure_law(law_name, law_params, nparams)
        CHARACTER(LEN=*), INTENT(IN) :: law_name
        INTEGER, INTENT(IN) :: nparams
        REAL*8, INTENT(IN), DIMENSION(nparams) :: law_params

        IF (.NOT. HOST_REGISTERED) THEN
            PRINT*, "WARNING: hostcluster not registered. Call add_hostcluster first"
            RETURN
        END IF

        IF (nparams < 1) THEN
            PRINT*, "WARNING: set_hostcluster_structure_law requires nparams >= 1"
            RETURN
        END IF

        IF (ALLOCATED(host_law_params)) DEALLOCATE(host_law_params)
        ALLOCATE(host_law_params(nparams))
        host_law_params = law_params

        HOST_LAW_NAME = TRIM(law_name)
        HOST_PARAM_MODE = HOST_PARAM_MODE_LAW
        HOST_FINALIZED = .FALSE.
    END SUBROUTINE set_hostcluster_structure_law

    SUBROUTINE set_hostcluster_structure_table(times, ptable, ntimes, nparams)
        INTEGER, INTENT(IN) :: ntimes, nparams
        REAL*8, INTENT(IN), DIMENSION(ntimes) :: times
        REAL*8, INTENT(IN), DIMENSION(ntimes, nparams) :: ptable

        IF (.NOT. HOST_REGISTERED) THEN
            PRINT*, "WARNING: hostcluster not registered. Call add_hostcluster first"
            RETURN
        END IF

        IF (ntimes < 2 .OR. nparams < 1) THEN
            PRINT*, "WARNING: set_hostcluster_structure_table requires ntimes>=2 and nparams>=1"
            RETURN
        END IF

        IF (ALLOCATED(host_param_times)) DEALLOCATE(host_param_times)
        IF (ALLOCATED(host_param_table)) DEALLOCATE(host_param_table)
        ALLOCATE(host_param_times(ntimes), host_param_table(ntimes, nparams))

        host_param_times = times
        host_param_table = ptable
        HOST_PARAM_MODE = HOST_PARAM_MODE_TABLE
        HOST_FINALIZED = .FALSE.
    END SUBROUTINE set_hostcluster_structure_table

    SUBROUTINE finalize_hostcluster()
        HOST_FINALIZED = .FALSE.

        IF (.NOT. HOST_REGISTERED) THEN
            PRINT*, "WARNING: finalize_hostcluster called before add_hostcluster"
            RETURN
        END IF

        IF (.NOT. HOST_KINEMATICS_SET) THEN
            PRINT*, "WARNING: finalize_hostcluster requires host kinematics"
            RETURN
        END IF

        IF (LEN_TRIM(HOST_BACKEND_NAME) < 1) THEN
            PRINT*, "WARNING: finalize_hostcluster requires backend selection"
            RETURN
        END IF

        IF (HOST_PARAM_MODE == HOST_PARAM_MODE_NONE) THEN
            PRINT*, "WARNING: finalize_hostcluster requires structure model"
            RETURN
        END IF

        HOST_FINALIZED = .TRUE.
    END SUBROUTINE finalize_hostcluster

    SUBROUTINE update_hostcluster_state(t)
        REAL*8, INTENT(IN) :: t

        IF (.NOT. HOST_FINALIZED) RETURN
        IF (.NOT. ALLOCATED(host_times)) RETURN

        CALL sample_kinematics_at_time(t)
    END SUBROUTINE update_hostcluster_state

    SUBROUTINE force_hostcluster_on_particles(nparticles, x, y, z, ax, ay, az, phi)
        INTEGER, INTENT(IN) :: nparticles
        REAL*8, INTENT(IN), DIMENSION(nparticles) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(nparticles) :: ax, ay, az, phi

        ! Base contract scaffold: backend physics is added in later phases.
        ax = 0.0D0
        ay = 0.0D0
        az = 0.0D0
        phi = 0.0D0
    END SUBROUTINE force_hostcluster_on_particles

    SUBROUTINE sample_kinematics_at_time(t)
        REAL*8, INTENT(IN) :: t
        INTEGER :: i, n
        REAL*8 :: alpha

        n = SIZE(host_times)
        IF (t <= host_times(1)) THEN
            host_x_current = host_x(1)
            host_y_current = host_y(1)
            host_z_current = host_z(1)
            host_vx_current = host_vx(1)
            host_vy_current = host_vy(1)
            host_vz_current = host_vz(1)
            RETURN
        END IF
        IF (t >= host_times(n)) THEN
            host_x_current = host_x(n)
            host_y_current = host_y(n)
            host_z_current = host_z(n)
            host_vx_current = host_vx(n)
            host_vy_current = host_vy(n)
            host_vz_current = host_vz(n)
            RETURN
        END IF

        DO i = 1, n - 1
            IF (host_times(i) <= t .AND. t <= host_times(i + 1)) THEN
                alpha = (t - host_times(i)) / (host_times(i + 1) - host_times(i))
                host_x_current = (1.0D0 - alpha) * host_x(i) + alpha * host_x(i + 1)
                host_y_current = (1.0D0 - alpha) * host_y(i) + alpha * host_y(i + 1)
                host_z_current = (1.0D0 - alpha) * host_z(i) + alpha * host_z(i + 1)
                host_vx_current = (1.0D0 - alpha) * host_vx(i) + alpha * host_vx(i + 1)
                host_vy_current = (1.0D0 - alpha) * host_vy(i) + alpha * host_vy(i + 1)
                host_vz_current = (1.0D0 - alpha) * host_vz(i) + alpha * host_vz(i + 1)
                RETURN
            END IF
        END DO
    END SUBROUTINE sample_kinematics_at_time

END MODULE hostcluster
