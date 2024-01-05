module constants
    ! set the gravitational constant and time step
    !-----------------------------------------------------------------------
    ! ALSO save the global time steps, which are described in the jupyter notebook TIMESTEP.IPYNB
    IMPLICIT NONE
    REAL*8,parameter :: G=4.300917270036279e-06 !! in solar masses and km/s
    REAL*8,parameter :: pi=3.141592653589793238462643 
    INTEGER,parameter:: saveskip=8 ! save every 8th step. savetimestep/minimumtimestep
    REAL*8,parameter :: minimumtimestep=2.5567804126142377e-06      ! in s kpc / km
    REAL*8,parameter :: savetimestep=2.04542433009139e-05           ! in s kpc / km
    REAL*8,parameter :: backwardintegrationtime=5.113560825228475   ! in s kpc / km
    REAL*8,parameter :: forwardintegrationtime=1.022712165045695    ! in s kpc / km


    contains

END MODULE
