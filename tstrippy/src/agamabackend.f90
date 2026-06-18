MODULE agamabackend
    IMPLICIT NONE

    INTEGER, PARAMETER :: MAX_AGAMA_COMPONENTS = 16

    CHARACTER(LEN=8), DIMENSION(MAX_AGAMA_COMPONENTS), SAVE :: AGAMA_HANDLES = ''
    CHARACTER(LEN=256), DIMENSION(MAX_AGAMA_COMPONENTS), SAVE :: AGAMA_SPECS = ''
    LOGICAL, DIMENSION(MAX_AGAMA_COMPONENTS), SAVE :: AGAMA_SPEC_IS_FILE = .FALSE.
    INTEGER, SAVE :: N_AGAMA_COMPONENTS = 0
    LOGICAL, SAVE :: AGAMA_FINALIZED = .FALSE.

    INTERFACE
        SUBROUTINE agama_initfromfile(c_obj, inifilename)
            CHARACTER(LEN=8), INTENT(OUT) :: c_obj
            CHARACTER(LEN=*), INTENT(IN) :: inifilename
        END SUBROUTINE agama_initfromfile

        SUBROUTINE agama_initfromparam(c_obj, params)
            CHARACTER(LEN=8), INTENT(OUT) :: c_obj
            CHARACTER(LEN=*), INTENT(IN) :: params
        END SUBROUTINE agama_initfromparam

        REAL*8 FUNCTION agama_potential(c_obj, xyz)
            CHARACTER(LEN=8), INTENT(IN) :: c_obj
            REAL*8, INTENT(IN), DIMENSION(3) :: xyz
        END FUNCTION agama_potential

        REAL*8 FUNCTION agama_potforce(c_obj, xyz, force)
            CHARACTER(LEN=8), INTENT(IN) :: c_obj
            REAL*8, INTENT(IN), DIMENSION(3) :: xyz
            REAL*8, INTENT(OUT), DIMENSION(3) :: force
        END FUNCTION agama_potforce

        REAL*8 FUNCTION agama_density(c_obj, xyz)
            CHARACTER(LEN=8), INTENT(IN) :: c_obj
            REAL*8, INTENT(IN), DIMENSION(3) :: xyz
        END FUNCTION agama_density
    END INTERFACE

CONTAINS

    FUNCTION agama_real_to_string(value) RESULT(text)
        REAL*8, INTENT(IN) :: value
        CHARACTER(LEN=64) :: text

        WRITE(text,'(G0)') value
    END FUNCTION agama_real_to_string

    SUBROUTINE clear()
        N_AGAMA_COMPONENTS = 0
        AGAMA_FINALIZED = .FALSE.
        AGAMA_HANDLES = ''
        AGAMA_SPECS = ''
        AGAMA_SPEC_IS_FILE = .FALSE.
    END SUBROUTINE clear

    SUBROUTINE add_component_from_param(params)
        CHARACTER(LEN=*), INTENT(IN) :: params

        IF (AGAMA_FINALIZED) THEN
            WRITE(*,'(A)') 'WARNING: add_component_from_param: cannot add components after finalize'
            RETURN
        END IF
        IF (N_AGAMA_COMPONENTS >= MAX_AGAMA_COMPONENTS) THEN
            WRITE(*,'(A)') 'WARNING: add_component_from_param: maximum number of Agama components reached'
            RETURN
        END IF

        N_AGAMA_COMPONENTS = N_AGAMA_COMPONENTS + 1
        AGAMA_SPECS(N_AGAMA_COMPONENTS) = ADJUSTL(params)
        AGAMA_SPEC_IS_FILE(N_AGAMA_COMPONENTS) = .FALSE.
    END SUBROUTINE add_component_from_param

    SUBROUTINE add_component_from_file(inifilename)
        CHARACTER(LEN=*), INTENT(IN) :: inifilename

        IF (AGAMA_FINALIZED) THEN
            WRITE(*,'(A)') 'WARNING: add_component_from_file: cannot add components after finalize'
            RETURN
        END IF
        IF (N_AGAMA_COMPONENTS >= MAX_AGAMA_COMPONENTS) THEN
            WRITE(*,'(A)') 'WARNING: add_component_from_file: maximum number of Agama components reached'
            RETURN
        END IF

        N_AGAMA_COMPONENTS = N_AGAMA_COMPONENTS + 1
        AGAMA_SPECS(N_AGAMA_COMPONENTS) = ADJUSTL(inifilename)
        AGAMA_SPEC_IS_FILE(N_AGAMA_COMPONENTS) = .TRUE.
    END SUBROUTINE add_component_from_file

    SUBROUTINE add_disk_component(sigma0, hR, hZ)
        REAL*8, INTENT(IN) :: sigma0, hR, hZ
        CHARACTER(LEN=256) :: spec

        spec = 'type=Disk surfaceDensity=' // TRIM(ADJUSTL(agama_real_to_string(sigma0))) // &
               ' scaleRadius=' // TRIM(ADJUSTL(agama_real_to_string(hR))) // &
               ' scaleHeight=' // TRIM(ADJUSTL(agama_real_to_string(hZ)))
        CALL add_component_from_param(spec)
    END SUBROUTINE add_disk_component

    SUBROUTINE finalize()
        INTEGER :: i

        IF (AGAMA_FINALIZED) RETURN

        DO i = 1, N_AGAMA_COMPONENTS
            IF (AGAMA_SPEC_IS_FILE(i)) THEN
                CALL agama_initfromfile(AGAMA_HANDLES(i), TRIM(AGAMA_SPECS(i)))
            ELSE
                CALL agama_initfromparam(AGAMA_HANDLES(i), TRIM(AGAMA_SPECS(i)))
            END IF
        END DO

        AGAMA_FINALIZED = .TRUE.
    END SUBROUTINE finalize

    SUBROUTINE ensure_finalized()
        IF (.NOT. AGAMA_FINALIZED) CALL finalize()
    END SUBROUTINE ensure_finalized

    INTEGER FUNCTION ncomponents()
        ncomponents = N_AGAMA_COMPONENTS
    END FUNCTION ncomponents

    INTEGER FUNCTION component_slot(i_comp)
        INTEGER, INTENT(IN) :: i_comp

        IF (i_comp < 1 .OR. i_comp > N_AGAMA_COMPONENTS) THEN
            component_slot = 0
        ELSE
            component_slot = i_comp
        END IF
    END FUNCTION component_slot

    SUBROUTINE force_component(i_comp, n, x, y, z, ax, ay, az)
        INTEGER, INTENT(IN) :: i_comp, n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: ax, ay, az
        INTEGER :: j, slot
        REAL*8, DIMENSION(3) :: xyz, force_tmp
        REAL*8 :: phi_tmp

        CALL ensure_finalized()
        ax = 0.0D0
        ay = 0.0D0
        az = 0.0D0
        slot = component_slot(i_comp)
        IF (slot == 0) RETURN

        DO j = 1, n
            xyz(1) = x(j)
            xyz(2) = y(j)
            xyz(3) = z(j)
            phi_tmp = agama_potforce(AGAMA_HANDLES(slot), xyz, force_tmp)
            ax(j) = force_tmp(1)
            ay(j) = force_tmp(2)
            az(j) = force_tmp(3)
        END DO
    END SUBROUTINE force_component

    SUBROUTINE potential_component(i_comp, n, x, y, z, phi)
        INTEGER, INTENT(IN) :: i_comp, n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: phi
        INTEGER :: j, slot
        REAL*8, DIMENSION(3) :: xyz

        CALL ensure_finalized()
        phi = 0.0D0
        slot = component_slot(i_comp)
        IF (slot == 0) RETURN

        DO j = 1, n
            xyz(1) = x(j)
            xyz(2) = y(j)
            xyz(3) = z(j)
            phi(j) = agama_potential(AGAMA_HANDLES(slot), xyz)
        END DO
    END SUBROUTINE potential_component

    SUBROUTINE density_component(i_comp, n, x, y, z, rho)
        INTEGER, INTENT(IN) :: i_comp, n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: rho
        INTEGER :: j, slot
        REAL*8, DIMENSION(3) :: xyz

        CALL ensure_finalized()
        rho = 0.0D0
        slot = component_slot(i_comp)
        IF (slot == 0) RETURN

        DO j = 1, n
            xyz(1) = x(j)
            xyz(2) = y(j)
            xyz(3) = z(j)
            rho(j) = agama_density(AGAMA_HANDLES(slot), xyz)
        END DO
    END SUBROUTINE density_component

    SUBROUTINE force(n, x, y, z, ax, ay, az)
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: ax, ay, az
        INTEGER :: i
        REAL*8, DIMENSION(n) :: ax_c, ay_c, az_c

        CALL ensure_finalized()
        ax = 0.0D0
        ay = 0.0D0
        az = 0.0D0

        DO i = 1, N_AGAMA_COMPONENTS
            CALL force_component(i, n, x, y, z, ax_c, ay_c, az_c)
            ax = ax + ax_c
            ay = ay + ay_c
            az = az + az_c
        END DO
    END SUBROUTINE force

    SUBROUTINE potential(n, x, y, z, phi)
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: phi
        INTEGER :: i
        REAL*8, DIMENSION(n) :: phi_c

        CALL ensure_finalized()
        phi = 0.0D0

        DO i = 1, N_AGAMA_COMPONENTS
            CALL potential_component(i, n, x, y, z, phi_c)
            phi = phi + phi_c
        END DO
    END SUBROUTINE potential

    SUBROUTINE density(n, x, y, z, rho)
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(n) :: rho
        INTEGER :: i
        REAL*8, DIMENSION(n) :: rho_c

        CALL ensure_finalized()
        rho = 0.0D0

        DO i = 1, N_AGAMA_COMPONENTS
            CALL density_component(i, n, x, y, z, rho_c)
            rho = rho + rho_c
        END DO
    END SUBROUTINE density

    SUBROUTINE agama_smoke_exponential_disk(phi, force_vec, ok)
        REAL*8, INTENT(OUT) :: phi
        REAL*8, INTENT(OUT), DIMENSION(3) :: force_vec
        LOGICAL, INTENT(OUT) :: ok

        REAL*8, DIMENSION(1) :: x, y, z, phi_arr, ax, ay, az

        CALL clear()
        CALL add_component_from_param('type=Disk surfaceDensity=1 scaleRadius=3 scaleHeight=0.3')
        CALL finalize()

        x(1) = 8.0D0
        y(1) = 0.0D0
        z(1) = 0.1D0

        CALL potential(1, x, y, z, phi_arr)
        CALL force(1, x, y, z, ax, ay, az)

        phi = phi_arr(1)
        force_vec(1) = ax(1)
        force_vec(2) = ay(1)
        force_vec(3) = az(1)
        ok = .TRUE.
    END SUBROUTINE agama_smoke_exponential_disk

    SUBROUTINE agama_smoke_from_file(inifilename, phi, force, ok)
        CHARACTER(LEN=*), INTENT(IN) :: inifilename
        REAL*8, INTENT(OUT) :: phi
        REAL*8, INTENT(OUT), DIMENSION(3) :: force
        LOGICAL, INTENT(OUT) :: ok

        REAL*8, DIMENSION(3) :: xyz

        CALL clear()
        CALL add_component_from_file(inifilename)
        CALL finalize()

        xyz(1) = 8.0D0
        xyz(2) = 0.0D0
        xyz(3) = 0.1D0

        phi = agama_potential(AGAMA_HANDLES(1), xyz)
        phi = agama_potforce(AGAMA_HANDLES(1), xyz, force)
        ok = .TRUE.
    END SUBROUTINE agama_smoke_from_file

END MODULE agamabackend
