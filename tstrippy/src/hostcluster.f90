MODULE hostcluster
    USE gravity, ONLY : plummerforce, plummerpotential
    USE mathutils, ONLY : linear_interp_scalar

    IMPLICIT NONE

    INTEGER, PARAMETER, PUBLIC :: HOST_BACKEND_NONE = 0
    INTEGER, PARAMETER, PUBLIC :: HOST_BACKEND_PLUMMER = 1
    INTEGER, PARAMETER, PUBLIC :: HOST_BACKEND_KING = 2

    INTEGER, PARAMETER, PUBLIC :: HOST_PARAM_MODEL_CONSTANT = 0
    INTEGER, PARAMETER, PUBLIC :: HOST_PARAM_MODEL_TABLE = 1

    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: xhost, yhost, zhost
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: vxhost, vyhost, vzhost
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: timehost

    REAL*8, PUBLIC :: xhostcurrent = 0.0D0
    REAL*8, PUBLIC :: yhostcurrent = 0.0D0
    REAL*8, PUBLIC :: zhostcurrent = 0.0D0
    REAL*8, PUBLIC :: vxhostcurrent = 0.0D0
    REAL*8, PUBLIC :: vyhostcurrent = 0.0D0
    REAL*8, PUBLIC :: vzhostcurrent = 0.0D0

    LOGICAL, PUBLIC :: host_times_decreasing = .FALSE.
    INTEGER, PUBLIC :: hosttimeindex = 1

    INTEGER, PUBLIC :: host_backend = HOST_BACKEND_NONE
    INTEGER, PUBLIC :: host_param_model = HOST_PARAM_MODEL_CONSTANT

    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_params_current
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_params_constant
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: host_param_times
    REAL*8, DIMENSION(:,:), PUBLIC, ALLOCATABLE :: host_param_table
    INTEGER, PUBLIC :: host_param_time_index = 1

    REAL*8, PUBLIC :: masshostcurrent = 0.0D0
    REAL*8, PUBLIC :: radiushostcurrent = 0.0D0

    LOGICAL, PUBLIC :: host_enabled = .FALSE.

    LOGICAL, PUBLIC :: tracking_enabled = .FALSE.
    REAL*8, PUBLIC :: binding_energy_epsilon = 0.0D0
    INTEGER, PUBLIC :: escape_confirm_steps = 1
    LOGICAL, DIMENSION(:), PUBLIC, ALLOCATABLE :: is_bound
    LOGICAL, DIMENSION(:), PUBLIC, ALLOCATABLE :: ever_escaped
    INTEGER, DIMENSION(:), PUBLIC, ALLOCATABLE :: consecutive_unbound
    REAL*8, DIMENSION(:), PUBLIC, ALLOCATABLE :: escape_time

    PUBLIC :: clear
    PUBLIC :: init_kinematics
    PUBLIC :: set_backend
    PUBLIC :: set_constant_params
    PUBLIC :: set_param_history
    PUBLIC :: update_state
    PUBLIC :: force_on_particles
    PUBLIC :: classify_bound
    PUBLIC :: configure_escape_tracking
    PUBLIC :: reset_escape_tracking

    ! Legacy compatibility names used by older call sites.
    PUBLIC :: host_init_kinematics, host_init_mass, host_init_radius
    PUBLIC :: findhosttimeindex, updatehoststate
    PUBLIC :: hostallocation, hostdeallocation, computeforcebyhosts

CONTAINS

    SUBROUTINE clear()
        CALL hostdeallocation()
        host_backend = HOST_BACKEND_NONE
        host_param_model = HOST_PARAM_MODEL_CONSTANT
        host_enabled = .FALSE.
        tracking_enabled = .FALSE.
        binding_energy_epsilon = 0.0D0
        escape_confirm_steps = 1
        xhostcurrent = 0.0D0
        yhostcurrent = 0.0D0
        zhostcurrent = 0.0D0
        vxhostcurrent = 0.0D0
        vyhostcurrent = 0.0D0
        vzhostcurrent = 0.0D0
        masshostcurrent = 0.0D0
        radiushostcurrent = 0.0D0
        hosttimeindex = 1
        host_param_time_index = 1
        host_times_decreasing = .FALSE.
    END SUBROUTINE clear

    SUBROUTINE init_kinematics(ntimesteps, t, x, y, z, vx, vy, vz)
        INTEGER, INTENT(IN) :: ntimesteps
        REAL*8, INTENT(IN), DIMENSION(ntimesteps) :: t, x, y, z, vx, vy, vz
        INTEGER :: i

        IF (ntimesteps < 2) THEN
            WRITE(*,'(A)') "WARNING: hostcluster.init_kinematics requires at least 2 time points"
            RETURN
        END IF

        IF (t(2) > t(1)) THEN
            host_times_decreasing = .FALSE.
            DO i = 2, ntimesteps
                IF (t(i) <= t(i-1)) THEN
                    WRITE(*,'(A)') "WARNING: hostcluster.init_kinematics requires strictly increasing or decreasing times"
                    RETURN
                END IF
            END DO
        ELSE IF (t(2) < t(1)) THEN
            host_times_decreasing = .TRUE.
            DO i = 2, ntimesteps
                IF (t(i) >= t(i-1)) THEN
                    WRITE(*,'(A)') "WARNING: hostcluster.init_kinematics requires strictly increasing or decreasing times"
                    RETURN
                END IF
            END DO
        ELSE
            WRITE(*,'(A)') "WARNING: hostcluster.init_kinematics requires strictly monotonic times"
            RETURN
        END IF

        CALL hostallocation(ntimesteps)
        xhost = x
        yhost = y
        zhost = z
        vxhost = vx
        vyhost = vy
        vzhost = vz
        timehost = t
        hosttimeindex = 1
        host_enabled = .TRUE.
    END SUBROUTINE init_kinematics

    SUBROUTINE set_backend(model_name)
        CHARACTER(LEN=*), INTENT(IN) :: model_name

        SELECT CASE (TRIM(model_name))
            CASE ("plummer")
                host_backend = HOST_BACKEND_PLUMMER
            CASE ("king")
                host_backend = HOST_BACKEND_KING
                WRITE(*,'(A)') "WARNING: hostcluster.set_backend('king') is registered but not implemented yet"
            CASE DEFAULT
                WRITE(*,'(A,A)') "WARNING: hostcluster.set_backend unknown model: ", TRIM(model_name)
                RETURN
        END SELECT
    END SUBROUTINE set_backend

    SUBROUTINE set_constant_params(params, nparams)
        INTEGER, INTENT(IN) :: nparams
        REAL*8, INTENT(IN), DIMENSION(nparams) :: params

        IF (nparams < 1) THEN
            WRITE(*,'(A)') "WARNING: hostcluster.set_constant_params requires at least 1 parameter"
            RETURN
        END IF

        IF (ALLOCATED(host_params_constant)) DEALLOCATE(host_params_constant)
        IF (ALLOCATED(host_params_current)) DEALLOCATE(host_params_current)

        ALLOCATE(host_params_constant(nparams))
        ALLOCATE(host_params_current(nparams))
        host_params_constant = params
        host_params_current = params
        host_param_model = HOST_PARAM_MODEL_CONSTANT

        IF (nparams >= 1) masshostcurrent = params(1)
        IF (nparams >= 2) radiushostcurrent = params(2)
    END SUBROUTINE set_constant_params

    SUBROUTINE set_param_history(ntimes, nparams, t, ptable)
        INTEGER, INTENT(IN) :: ntimes, nparams
        REAL*8, INTENT(IN), DIMENSION(ntimes) :: t
        REAL*8, INTENT(IN), DIMENSION(ntimes, nparams) :: ptable
        INTEGER :: i

        IF (ntimes < 2) THEN
            WRITE(*,'(A)') "WARNING: hostcluster.set_param_history requires at least 2 time points"
            RETURN
        END IF
        IF (nparams < 1) THEN
            WRITE(*,'(A)') "WARNING: hostcluster.set_param_history requires at least 1 parameter"
            RETURN
        END IF

        IF (t(2) > t(1)) THEN
            DO i = 2, ntimes
                IF (t(i) <= t(i-1)) THEN
                    WRITE(*,'(A)') "WARNING: hostcluster.set_param_history times must be strictly monotonic"
                    RETURN
                END IF
            END DO
        ELSE IF (t(2) < t(1)) THEN
            DO i = 2, ntimes
                IF (t(i) >= t(i-1)) THEN
                    WRITE(*,'(A)') "WARNING: hostcluster.set_param_history times must be strictly monotonic"
                    RETURN
                END IF
            END DO
        ELSE
            WRITE(*,'(A)') "WARNING: hostcluster.set_param_history times must be strictly monotonic"
            RETURN
        END IF

        IF (ALLOCATED(host_param_times)) DEALLOCATE(host_param_times)
        IF (ALLOCATED(host_param_table)) DEALLOCATE(host_param_table)
        IF (ALLOCATED(host_params_current)) DEALLOCATE(host_params_current)

        ALLOCATE(host_param_times(ntimes))
        ALLOCATE(host_param_table(ntimes, nparams))
        ALLOCATE(host_params_current(nparams))

        host_param_times = t
        host_param_table = ptable
        host_params_current = ptable(1,:)
        host_param_model = HOST_PARAM_MODEL_TABLE
        host_param_time_index = 1

        masshostcurrent = host_params_current(1)
        IF (nparams >= 2) radiushostcurrent = host_params_current(2)
    END SUBROUTINE set_param_history

    SUBROUTINE update_state(mytime)
        REAL*8, INTENT(IN) :: mytime
        REAL*8 :: alpha

        IF (.NOT. ALLOCATED(timehost)) THEN
            WRITE(*,'(A)') "WARNING: hostcluster.update_state called before host kinematics were initialized"
            RETURN
        END IF

        CALL compute_interp_alpha(timehost, host_times_decreasing, hosttimeindex, mytime, alpha)

        xhostcurrent = linear_interp_scalar(xhost(hosttimeindex), xhost(hosttimeindex + 1), alpha)
        yhostcurrent = linear_interp_scalar(yhost(hosttimeindex), yhost(hosttimeindex + 1), alpha)
        zhostcurrent = linear_interp_scalar(zhost(hosttimeindex), zhost(hosttimeindex + 1), alpha)
        vxhostcurrent = linear_interp_scalar(vxhost(hosttimeindex), vxhost(hosttimeindex + 1), alpha)
        vyhostcurrent = linear_interp_scalar(vyhost(hosttimeindex), vyhost(hosttimeindex + 1), alpha)
        vzhostcurrent = linear_interp_scalar(vzhost(hosttimeindex), vzhost(hosttimeindex + 1), alpha)

        CALL update_current_params(mytime)
    END SUBROUTINE update_state

    SUBROUTINE force_on_particles(nparticles, x, y, z, ax, ay, az, phi)
        INTEGER, INTENT(IN) :: nparticles
        REAL*8, INTENT(IN), DIMENSION(nparticles) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(nparticles) :: ax, ay, az, phi

        REAL*8, DIMENSION(nparticles) :: dx, dy, dz, phihost
        REAL*8, DIMENSION(nparticles,3) :: forcehost
        REAL*8, DIMENSION(2) :: params

        ax = 0.0D0
        ay = 0.0D0
        az = 0.0D0
        phi = 0.0D0

        IF (.NOT. host_enabled) RETURN

        SELECT CASE (host_backend)
            CASE (HOST_BACKEND_PLUMMER)
                IF (.NOT. ALLOCATED(host_params_current)) THEN
                    WRITE(*,'(A)') "WARNING: hostcluster.force_on_particles has no host parameters"
                    RETURN
                END IF
                IF (SIZE(host_params_current) < 2) THEN
                    WRITE(*,'(A)') "WARNING: hostcluster.force_on_particles(plummer) requires params [M, a]"
                    RETURN
                END IF

                dx = x - xhostcurrent
                dy = y - yhostcurrent
                dz = z - zhostcurrent
                params(1) = host_params_current(1)
                params(2) = host_params_current(2)

                CALL plummerforce(params, nparticles, dx, dy, dz, forcehost)
                CALL plummerpotential(params, nparticles, dx, dy, dz, phihost)

                ax = forcehost(:,1)
                ay = forcehost(:,2)
                az = forcehost(:,3)
                phi = phihost

            CASE (HOST_BACKEND_KING)
                WRITE(*,'(A)') "WARNING: hostcluster.force_on_particles king backend is not implemented yet"
            CASE DEFAULT
                WRITE(*,'(A)') "WARNING: hostcluster.force_on_particles called before backend was configured"
        END SELECT
    END SUBROUTINE force_on_particles

    SUBROUTINE classify_bound(nparticles, x, y, z, vx, vy, vz, bound_out, erel_out)
        INTEGER, INTENT(IN) :: nparticles
        REAL*8, INTENT(IN), DIMENSION(nparticles) :: x, y, z, vx, vy, vz
        LOGICAL, INTENT(OUT), DIMENSION(nparticles) :: bound_out
        REAL*8, INTENT(OUT), DIMENSION(nparticles) :: erel_out

        REAL*8, DIMENSION(nparticles) :: dx, dy, dz, dvx, dvy, dvz
        REAL*8, DIMENSION(nparticles) :: ax_dummy, ay_dummy, az_dummy, phihost
        REAL*8, DIMENSION(nparticles) :: v2

        bound_out = .FALSE.
        erel_out = 0.0D0

        IF (.NOT. host_enabled) RETURN

        dx = x - xhostcurrent
        dy = y - yhostcurrent
        dz = z - zhostcurrent
        dvx = vx - vxhostcurrent
        dvy = vy - vyhostcurrent
        dvz = vz - vzhostcurrent
        v2 = dvx*dvx + dvy*dvy + dvz*dvz

        CALL force_on_particles(nparticles, x, y, z, ax_dummy, ay_dummy, az_dummy, phihost)
        erel_out = 0.5D0 * v2 + phihost
        bound_out = erel_out < (-binding_energy_epsilon)
    END SUBROUTINE classify_bound

    SUBROUTINE configure_escape_tracking(nparticles, epsilon, confirm_steps)
        INTEGER, INTENT(IN) :: nparticles
        REAL*8, INTENT(IN) :: epsilon
        INTEGER, INTENT(IN) :: confirm_steps

        IF (nparticles < 1) THEN
            WRITE(*,'(A)') "WARNING: hostcluster.configure_escape_tracking nparticles must be >= 1"
            RETURN
        END IF
        IF (confirm_steps < 1) THEN
            WRITE(*,'(A)') "WARNING: hostcluster.configure_escape_tracking confirm_steps must be >= 1"
            RETURN
        END IF

        IF (ALLOCATED(is_bound)) DEALLOCATE(is_bound)
        IF (ALLOCATED(ever_escaped)) DEALLOCATE(ever_escaped)
        IF (ALLOCATED(consecutive_unbound)) DEALLOCATE(consecutive_unbound)
        IF (ALLOCATED(escape_time)) DEALLOCATE(escape_time)

        ALLOCATE(is_bound(nparticles))
        ALLOCATE(ever_escaped(nparticles))
        ALLOCATE(consecutive_unbound(nparticles))
        ALLOCATE(escape_time(nparticles))

        is_bound = .TRUE.
        ever_escaped = .FALSE.
        consecutive_unbound = 0
        escape_time = -1.0D0

        binding_energy_epsilon = epsilon
        escape_confirm_steps = confirm_steps
        tracking_enabled = .TRUE.
    END SUBROUTINE configure_escape_tracking

    SUBROUTINE reset_escape_tracking()
        IF (ALLOCATED(is_bound)) is_bound = .TRUE.
        IF (ALLOCATED(ever_escaped)) ever_escaped = .FALSE.
        IF (ALLOCATED(consecutive_unbound)) consecutive_unbound = 0
        IF (ALLOCATED(escape_time)) escape_time = -1.0D0
    END SUBROUTINE reset_escape_tracking

    SUBROUTINE update_current_params(mytime)
        REAL*8, INTENT(IN) :: mytime
        REAL*8 :: alpha

        IF (.NOT. ALLOCATED(host_params_current)) RETURN

        SELECT CASE (host_param_model)
            CASE (HOST_PARAM_MODEL_CONSTANT)
                IF (ALLOCATED(host_params_constant)) host_params_current = host_params_constant
            CASE (HOST_PARAM_MODEL_TABLE)
                IF (.NOT. ALLOCATED(host_param_times)) RETURN
                CALL compute_interp_alpha(host_param_times, param_times_decreasing(), host_param_time_index, mytime, alpha)
                host_params_current = linear_interp_vector(
                    host_param_table(host_param_time_index,:), &
                    host_param_table(host_param_time_index + 1,:), alpha)
            CASE DEFAULT
                RETURN
        END SELECT

        IF (SIZE(host_params_current) >= 1) masshostcurrent = host_params_current(1)
        IF (SIZE(host_params_current) >= 2) radiushostcurrent = host_params_current(2)
    END SUBROUTINE update_current_params

    LOGICAL FUNCTION param_times_decreasing()
        IF (.NOT. ALLOCATED(host_param_times)) THEN
            param_times_decreasing = .FALSE.
            RETURN
        END IF
        param_times_decreasing = (host_param_times(2) < host_param_times(1))
    END FUNCTION param_times_decreasing

    SUBROUTINE compute_interp_alpha(tarr, decreasing, idx, tvalue, alpha)
        REAL*8, INTENT(IN), DIMENSION(:) :: tarr
        LOGICAL, INTENT(IN) :: decreasing
        INTEGER, INTENT(INOUT) :: idx
        REAL*8, INTENT(IN) :: tvalue
        REAL*8, INTENT(OUT) :: alpha

        INTEGER :: ntime
        REAL*8 :: t0, t1, denom

        ntime = SIZE(tarr)
        IF (ntime < 2) THEN
            alpha = 0.0D0
            RETURN
        END IF

        IF (.NOT. decreasing) THEN
            IF (tvalue <= tarr(1)) THEN
                idx = 1
                alpha = 0.0D0
                RETURN
            ELSE IF (tvalue >= tarr(ntime)) THEN
                idx = ntime - 1
                alpha = 1.0D0
                RETURN
            END IF

            DO WHILE (idx < ntime - 1 .AND. tarr(idx + 1) < tvalue)
                idx = idx + 1
            END DO
            DO WHILE (idx > 1 .AND. tarr(idx) > tvalue)
                idx = idx - 1
            END DO
        ELSE
            IF (tvalue >= tarr(1)) THEN
                idx = 1
                alpha = 0.0D0
                RETURN
            ELSE IF (tvalue <= tarr(ntime)) THEN
                idx = ntime - 1
                alpha = 1.0D0
                RETURN
            END IF

            DO WHILE (idx < ntime - 1 .AND. tarr(idx + 1) > tvalue)
                idx = idx + 1
            END DO
            DO WHILE (idx > 1 .AND. tarr(idx) < tvalue)
                idx = idx - 1
            END DO
        END IF

        t0 = tarr(idx)
        t1 = tarr(idx + 1)
        denom = t1 - t0
        IF (denom == 0.0D0) THEN
            alpha = 0.0D0
        ELSE
            alpha = (tvalue - t0) / denom
        END IF
    END SUBROUTINE compute_interp_alpha

    FUNCTION linear_interp_vector(v0, v1, alpha) RESULT(vout)
        REAL*8, INTENT(IN), DIMENSION(:) :: v0, v1
        REAL*8, INTENT(IN) :: alpha
        REAL*8, DIMENSION(SIZE(v0)) :: vout

        vout = (1.0D0 - alpha) * v0 + alpha * v1
    END FUNCTION linear_interp_vector

    SUBROUTINE host_init_kinematics(NTIMESTEPS, t, x, y, z, vx, vy, vz)
        INTEGER, INTENT(IN) :: NTIMESTEPS
        REAL*8, INTENT(IN), DIMENSION(NTIMESTEPS) :: t, x, y, z, vx, vy, vz
        CALL init_kinematics(NTIMESTEPS, t, x, y, z, vx, vy, vz)
    END SUBROUTINE host_init_kinematics

    SUBROUTINE host_init_mass(mass_model, params)
        INTEGER, INTENT(IN) :: mass_model
        REAL*8, INTENT(IN), DIMENSION(:) :: params

        SELECT CASE (mass_model)
            CASE (0)
                CALL set_constant_params(params, SIZE(params))
            CASE DEFAULT
                WRITE(*,'(A)') "WARNING: hostcluster.host_init_mass only supports constant model in Phase 1"
                CALL set_constant_params(params, SIZE(params))
        END SELECT
    END SUBROUTINE host_init_mass

    SUBROUTINE host_init_radius(radius)
        REAL*8, INTENT(IN) :: radius
        IF (.NOT. ALLOCATED(host_params_current)) THEN
            ALLOCATE(host_params_current(2))
            host_params_current = 0.0D0
        END IF
        IF (SIZE(host_params_current) < 2) THEN
            WRITE(*,'(A)') "WARNING: hostcluster.host_init_radius expects at least 2 host params [M, a]"
            RETURN
        END IF
        host_params_current(2) = radius
        radiushostcurrent = radius
    END SUBROUTINE host_init_radius

    SUBROUTINE updatehoststate(mytime)
        REAL*8, INTENT(IN) :: mytime
        CALL update_state(mytime)
    END SUBROUTINE updatehoststate

    SUBROUTINE findhosttimeindex(mytime)
        REAL*8, INTENT(IN) :: mytime
        CALL update_state(mytime)
    END SUBROUTINE findhosttimeindex

    SUBROUTINE computeforcebyhosts(Nparticles, x, y, z, ax, ay, az, phi)
        INTEGER, INTENT(IN) :: Nparticles
        REAL*8, INTENT(IN), DIMENSION(Nparticles) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(Nparticles) :: ax, ay, az, phi
        CALL force_on_particles(Nparticles, x, y, z, ax, ay, az, phi)
    END SUBROUTINE computeforcebyhosts

    SUBROUTINE hostallocation(NTIMESTEPS)
        INTEGER, INTENT(IN) :: NTIMESTEPS
        IF (ALLOCATED(xhost)) DEALLOCATE(xhost)
        IF (ALLOCATED(yhost)) DEALLOCATE(yhost)
        IF (ALLOCATED(zhost)) DEALLOCATE(zhost)
        IF (ALLOCATED(vxhost)) DEALLOCATE(vxhost)
        IF (ALLOCATED(vyhost)) DEALLOCATE(vyhost)
        IF (ALLOCATED(vzhost)) DEALLOCATE(vzhost)
        IF (ALLOCATED(timehost)) DEALLOCATE(timehost)

        ALLOCATE(xhost(NTIMESTEPS))
        ALLOCATE(yhost(NTIMESTEPS))
        ALLOCATE(zhost(NTIMESTEPS))
        ALLOCATE(vxhost(NTIMESTEPS))
        ALLOCATE(vyhost(NTIMESTEPS))
        ALLOCATE(vzhost(NTIMESTEPS))
        ALLOCATE(timehost(NTIMESTEPS))
    END SUBROUTINE hostallocation

    SUBROUTINE hostdeallocation
        IF (ALLOCATED(xhost)) DEALLOCATE(xhost)
        IF (ALLOCATED(yhost)) DEALLOCATE(yhost)
        IF (ALLOCATED(zhost)) DEALLOCATE(zhost)
        IF (ALLOCATED(vxhost)) DEALLOCATE(vxhost)
        IF (ALLOCATED(vyhost)) DEALLOCATE(vyhost)
        IF (ALLOCATED(vzhost)) DEALLOCATE(vzhost)
        IF (ALLOCATED(timehost)) DEALLOCATE(timehost)
        IF (ALLOCATED(host_params_current)) DEALLOCATE(host_params_current)
        IF (ALLOCATED(host_params_constant)) DEALLOCATE(host_params_constant)
        IF (ALLOCATED(host_param_times)) DEALLOCATE(host_param_times)
        IF (ALLOCATED(host_param_table)) DEALLOCATE(host_param_table)
        IF (ALLOCATED(is_bound)) DEALLOCATE(is_bound)
        IF (ALLOCATED(ever_escaped)) DEALLOCATE(ever_escaped)
        IF (ALLOCATED(consecutive_unbound)) DEALLOCATE(consecutive_unbound)
        IF (ALLOCATED(escape_time)) DEALLOCATE(escape_time)
    END SUBROUTINE hostdeallocation

END MODULE hostcluster
