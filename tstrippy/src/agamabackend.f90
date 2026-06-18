MODULE agamabackend
    IMPLICIT NONE

CONTAINS

    SUBROUTINE agama_smoke_exponential_disk(phi, force, ok)
        ! Minimal smoke call that constructs one Agama Disk potential and evaluates force at one point.
        REAL*8, INTENT(OUT) :: phi
        REAL*8, INTENT(OUT), DIMENSION(3) :: force
        LOGICAL, INTENT(OUT) :: ok

        CHARACTER(LEN=8) :: c_obj
        REAL*8, DIMENSION(3) :: xyz
        REAL*8 :: agama_potforce
        EXTERNAL :: agama_initfromparam
        EXTERNAL :: agama_potforce

        phi = 0.0D0
        force = 0.0D0
        ok = .FALSE.

        CALL agama_initfromparam(c_obj, &
            'type=Disk surfaceDensity=1 scaleRadius=3 scaleHeight=0.3')

        xyz(1) = 8.0D0
        xyz(2) = 0.0D0
        xyz(3) = 0.1D0

        phi = agama_potforce(c_obj, xyz, force)
        ok = .TRUE.
    END SUBROUTINE agama_smoke_exponential_disk

END MODULE agamabackend
