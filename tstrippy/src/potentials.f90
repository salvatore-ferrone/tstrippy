MODULE potentials
    USE mathutils, ONLY: linear_interp_scalar, legendre_axisymmetric_basis, &
                         legendre_p_all_axisymmetric, gauss_legendre_nodes_weights, &
                         bessel_j0_scalar, bessel_j1_scalar
    IMPLICIT NONE

    ABSTRACT INTERFACE
        SUBROUTINE bessel_component_evaluator_interface(params, N, x, y, z, ax, ay, az, phi)
            REAL*8, INTENT(IN), DIMENSION(*) :: params
            INTEGER, INTENT(IN) :: N
            REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
            REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi
        END SUBROUTINE bessel_component_evaluator_interface
    END INTERFACE

    ! Axisymmetric basis-expansion tables (l-mode by radial grid index).
    ! BASIS_GRID_SET: r_grid allocated, arrays ready for projection.
    ! BASIS_EXPANSION_INITIALIZED: projection + potential tables computed, ready to evaluate.
    LOGICAL, PUBLIC  :: BASIS_GRID_SET              = .FALSE.
    LOGICAL, PUBLIC  :: BASIS_EXPANSION_INITIALIZED = .FALSE.
    INTEGER, PUBLIC  :: BASIS_LMAX = -1
    INTEGER, PUBLIC  :: BASIS_NR   = -1
    REAL*8,  PUBLIC  :: BASIS_G    = -1.0D0

    REAL*8, DIMENSION(:),   ALLOCATABLE, PUBLIC :: BASIS_R_GRID
    REAL*8, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: BASIS_RHO_L_GRID
    REAL*8, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: BASIS_PHI_L_GRID
    REAL*8, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: BASIS_DPHI_L_DR_GRID

    ! Composite BFE manager state.
    INTEGER, PARAMETER, PUBLIC :: BFE_KIND_NONE = 0
    INTEGER, PARAMETER, PUBLIC :: BFE_KIND_LEGENDRE = 1
    INTEGER, PARAMETER, PUBLIC :: BFE_KIND_BESSEL_DISK = 2
    INTEGER, PARAMETER, PUBLIC :: COMPOSITE_BESSEL_MAX_PARAMS = 16

    LOGICAL, PUBLIC :: COMPOSITE_BASIS_GRID_SET = .FALSE.
    LOGICAL, PUBLIC :: COMPOSITE_BASIS_FINALIZED = .FALSE.
    INTEGER, PUBLIC :: COMPOSITE_NCOMP = 0
    INTEGER, PUBLIC :: COMPOSITE_LMAX = -1
    INTEGER, PUBLIC :: COMPOSITE_NR = -1
    REAL*8, PUBLIC :: COMPOSITE_G = -1.0D0

    INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: COMPOSITE_KIND
    LOGICAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: COMPOSITE_READY
    INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: COMPOSITE_BESSEL_NPARAMS

    REAL*8, DIMENSION(:), ALLOCATABLE, PUBLIC :: COMPOSITE_R_GRID
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: COMPOSITE_RHO_L_GRID
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: COMPOSITE_PHI_L_GRID
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: COMPOSITE_DPHI_L_DR_GRID

    REAL*8, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: COMPOSITE_BESSEL_PARAMS

    ! Private internal routines (not exposed to Python/caller)
    PRIVATE :: default_init_basis_expansion
    PRIVATE :: project_exponential_oblate_halo
    PRIVATE :: compute_phi_tables_from_rho
    PRIVATE :: axisymmetricbasisexpansion_eval
    PRIVATE :: project_ibata2024halo
    PRIVATE :: compute_phi_tables_from_rho_component
    PRIVATE :: axisymmetricbasisexpansion_eval_component
    PRIVATE :: bessel_eval_component
    PRIVATE :: exponential_disk_bessel_eval_component

    CONTAINS
    SUBROUTINE initaxisymmetricbasisexpansion(G, lmax, nr, r_grid)
        ! Allocate basis-expansion storage and store the radial grid.
        ! Does NOT project any density or compute potential tables.
        ! Call a density subroutine (e.g. exponential_oblate_halo) afterwards
        ! to trigger projection and table construction.
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: G
        INTEGER, INTENT(IN) :: lmax, nr
        REAL*8, INTENT(IN), DIMENSION(nr) :: r_grid

        CALL clearaxisymmetricbasisexpansion()

        ALLOCATE(BASIS_R_GRID(nr))
        ALLOCATE(BASIS_RHO_L_GRID(0:lmax, nr))
        ALLOCATE(BASIS_PHI_L_GRID(0:lmax, nr))
        ALLOCATE(BASIS_DPHI_L_DR_GRID(0:lmax, nr))

        BASIS_RHO_L_GRID     = 0.0D0
        BASIS_PHI_L_GRID     = 0.0D0
        BASIS_DPHI_L_DR_GRID = 0.0D0

        BASIS_G    = G
        BASIS_LMAX = lmax
        BASIS_NR   = nr
        BASIS_R_GRID = r_grid

        BASIS_GRID_SET              = .TRUE.
        BASIS_EXPANSION_INITIALIZED = .FALSE.
    END SUBROUTINE initaxisymmetricbasisexpansion

    SUBROUTINE clearaxisymmetricbasisexpansion()
        IMPLICIT NONE

        IF (ALLOCATED(BASIS_R_GRID))          DEALLOCATE(BASIS_R_GRID)
        IF (ALLOCATED(BASIS_RHO_L_GRID))      DEALLOCATE(BASIS_RHO_L_GRID)
        IF (ALLOCATED(BASIS_PHI_L_GRID))      DEALLOCATE(BASIS_PHI_L_GRID)
        IF (ALLOCATED(BASIS_DPHI_L_DR_GRID))  DEALLOCATE(BASIS_DPHI_L_DR_GRID)

        BASIS_GRID_SET              = .FALSE.
        BASIS_EXPANSION_INITIALIZED = .FALSE.
        BASIS_LMAX = -1
        BASIS_NR   = -1
        BASIS_G    = -1.0D0
    END SUBROUTINE clearaxisymmetricbasisexpansion

    SUBROUTINE clearaxisymmetriccompositebasisexpansion()
        IMPLICIT NONE

        IF (ALLOCATED(COMPOSITE_KIND))            DEALLOCATE(COMPOSITE_KIND)
        IF (ALLOCATED(COMPOSITE_READY))           DEALLOCATE(COMPOSITE_READY)
        IF (ALLOCATED(COMPOSITE_R_GRID))          DEALLOCATE(COMPOSITE_R_GRID)
        IF (ALLOCATED(COMPOSITE_RHO_L_GRID))      DEALLOCATE(COMPOSITE_RHO_L_GRID)
        IF (ALLOCATED(COMPOSITE_PHI_L_GRID))      DEALLOCATE(COMPOSITE_PHI_L_GRID)
        IF (ALLOCATED(COMPOSITE_DPHI_L_DR_GRID))  DEALLOCATE(COMPOSITE_DPHI_L_DR_GRID)
        IF (ALLOCATED(COMPOSITE_BESSEL_NPARAMS))  DEALLOCATE(COMPOSITE_BESSEL_NPARAMS)
        IF (ALLOCATED(COMPOSITE_BESSEL_PARAMS))   DEALLOCATE(COMPOSITE_BESSEL_PARAMS)

        COMPOSITE_BASIS_GRID_SET = .FALSE.
        COMPOSITE_BASIS_FINALIZED = .FALSE.
        COMPOSITE_NCOMP = 0
        COMPOSITE_LMAX = -1
        COMPOSITE_NR = -1
        COMPOSITE_G = -1.0D0
    END SUBROUTINE clearaxisymmetriccompositebasisexpansion

    SUBROUTINE clearcompositebasisexpansion()
        ! Backward-compatible alias.
        IMPLICIT NONE
        CALL clearaxisymmetriccompositebasisexpansion()
    END SUBROUTINE clearcompositebasisexpansion

    SUBROUTINE initaxisymmetriccompositebasisexpansion(G, lmax, nr, r_grid, ncomp)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: G
        INTEGER, INTENT(IN) :: lmax, nr, ncomp
        REAL*8, INTENT(IN), DIMENSION(nr) :: r_grid

        CALL clearaxisymmetriccompositebasisexpansion()

        ALLOCATE(COMPOSITE_KIND(ncomp), COMPOSITE_READY(ncomp))
        ALLOCATE(COMPOSITE_BESSEL_NPARAMS(ncomp))
        ALLOCATE(COMPOSITE_R_GRID(nr))
        ALLOCATE(COMPOSITE_RHO_L_GRID(0:lmax, nr, ncomp))
        ALLOCATE(COMPOSITE_PHI_L_GRID(0:lmax, nr, ncomp))
        ALLOCATE(COMPOSITE_DPHI_L_DR_GRID(0:lmax, nr, ncomp))
        ALLOCATE(COMPOSITE_BESSEL_PARAMS(COMPOSITE_BESSEL_MAX_PARAMS, ncomp))

        COMPOSITE_KIND = BFE_KIND_NONE
        COMPOSITE_READY = .FALSE.
        COMPOSITE_BESSEL_NPARAMS = 0
        COMPOSITE_RHO_L_GRID = 0.0D0
        COMPOSITE_PHI_L_GRID = 0.0D0
        COMPOSITE_DPHI_L_DR_GRID = 0.0D0
        COMPOSITE_BESSEL_PARAMS = 0.0D0

        COMPOSITE_G = G
        COMPOSITE_LMAX = lmax
        COMPOSITE_NR = nr
        COMPOSITE_NCOMP = ncomp
        COMPOSITE_R_GRID = r_grid
        COMPOSITE_BASIS_GRID_SET = .TRUE.
        COMPOSITE_BASIS_FINALIZED = .FALSE.
    END SUBROUTINE initaxisymmetriccompositebasisexpansion

    SUBROUTINE initcompositebasisexpansion(G, lmax, nr, r_grid, ncomp)
        ! Backward-compatible alias.
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: G
        INTEGER, INTENT(IN) :: lmax, nr, ncomp
        REAL*8, INTENT(IN), DIMENSION(nr) :: r_grid

        CALL initaxisymmetriccompositebasisexpansion(G, lmax, nr, r_grid, ncomp)
    END SUBROUTINE initcompositebasisexpansion

    SUBROUTINE addcompositeexponentialoblate(component_index, rho0, s0, q)
        ! Legendre basis path: preferred for spherical or near-spherical components.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: component_index
        REAL*8, INTENT(IN) :: rho0, s0, q
        INTEGER :: n_mu, i_r, l, k
        REAL*8, ALLOCATABLE :: mu_q(:), w_q(:), p(:)
        REAL*8 :: mu, rho_val, eta, factor

        IF (.NOT. COMPOSITE_BASIS_GRID_SET) STOP "initaxisymmetriccompositebasisexpansion must be called before addcompositeexponentialoblate"
        IF (component_index < 1 .OR. component_index > COMPOSITE_NCOMP) STOP "invalid composite component_index"

        n_mu = MAX(4 * (COMPOSITE_LMAX + 1), 40)
        ALLOCATE(mu_q(n_mu), w_q(n_mu), p(0:COMPOSITE_LMAX))
        CALL gauss_legendre_nodes_weights(n_mu, mu_q, w_q)

        eta = 1.0D0 - 1.0D0 / (q * q)
        COMPOSITE_RHO_L_GRID(:,:,component_index) = 0.0D0
        DO i_r = 1, COMPOSITE_NR
            DO k = 1, n_mu
                mu = mu_q(k)
                rho_val = rho0 * EXP(-COMPOSITE_R_GRID(i_r) / s0 * SQRT(MAX(1.0D0 - eta*mu*mu, 0.0D0)))
                CALL legendre_p_all_axisymmetric(COMPOSITE_LMAX, mu, p)
                DO l = 0, COMPOSITE_LMAX, 2
                    factor = (2*l + 1) * 0.5D0 * w_q(k)
                    COMPOSITE_RHO_L_GRID(l, i_r, component_index) = COMPOSITE_RHO_L_GRID(l, i_r, component_index) + factor * rho_val * p(l)
                END DO
            END DO
        END DO
        DEALLOCATE(mu_q, w_q, p)

        CALL compute_phi_tables_from_rho_component(COMPOSITE_R_GRID, COMPOSITE_G, &
            COMPOSITE_RHO_L_GRID(:,:,component_index), COMPOSITE_PHI_L_GRID(:,:,component_index), &
            COMPOSITE_DPHI_L_DR_GRID(:,:,component_index))

        COMPOSITE_KIND(component_index) = BFE_KIND_LEGENDRE
        COMPOSITE_READY(component_index) = .TRUE.
        COMPOSITE_BASIS_FINALIZED = .FALSE.
    END SUBROUTINE addcompositeexponentialoblate

    SUBROUTINE addcompositeibata2024halo(component_index, rho0, r0, rt, q, gamma, beta)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: component_index
        REAL*8, INTENT(IN) :: rho0, r0, rt, q, gamma, beta
        INTEGER :: n_mu, i_r, l, k
        REAL*8, ALLOCATABLE :: mu_q(:), w_q(:), p(:)
        REAL*8 :: mu, rho_val, eta, factor, s, x
        REAL*8, PARAMETER :: s_floor = 1.0D-12

        IF (.NOT. COMPOSITE_BASIS_GRID_SET) STOP "initaxisymmetriccompositebasisexpansion must be called before addcompositeibata2024"
        IF (component_index < 1 .OR. component_index > COMPOSITE_NCOMP) STOP "invalid composite component_index"

        n_mu = MAX(4 * (COMPOSITE_LMAX + 1), 40)
        ALLOCATE(mu_q(n_mu), w_q(n_mu), p(0:COMPOSITE_LMAX))
        CALL gauss_legendre_nodes_weights(n_mu, mu_q, w_q)

        eta = 1.0D0 - 1.0D0 / (q * q)
        COMPOSITE_RHO_L_GRID(:,:,component_index) = 0.0D0
        DO i_r = 1, COMPOSITE_NR
            DO k = 1, n_mu
                mu = mu_q(k)
                s = COMPOSITE_R_GRID(i_r) * SQRT(MAX(1.0D0 - eta*mu*mu, 0.0D0))
                x = MAX(s / r0, s_floor)
                rho_val = rho0 * x**(-gamma) * (1.0D0 + x)**(gamma - beta) * EXP(-(s/rt)**2)
                CALL legendre_p_all_axisymmetric(COMPOSITE_LMAX, mu, p)
                DO l = 0, COMPOSITE_LMAX, 2
                    factor = (2*l + 1) * 0.5D0 * w_q(k)
                    COMPOSITE_RHO_L_GRID(l, i_r, component_index) = COMPOSITE_RHO_L_GRID(l, i_r, component_index) + factor * rho_val * p(l)
                END DO
            END DO
        END DO
        DEALLOCATE(mu_q, w_q, p)

        CALL compute_phi_tables_from_rho_component(COMPOSITE_R_GRID, COMPOSITE_G, &
            COMPOSITE_RHO_L_GRID(:,:,component_index), COMPOSITE_PHI_L_GRID(:,:,component_index), &
            COMPOSITE_DPHI_L_DR_GRID(:,:,component_index))

        COMPOSITE_KIND(component_index) = BFE_KIND_LEGENDRE
        COMPOSITE_READY(component_index) = .TRUE.
        COMPOSITE_BASIS_FINALIZED = .FALSE.
    END SUBROUTINE addcompositeibata2024halo

    SUBROUTINE addcompositebesselcomponent(component_index, params, nparams)
        ! Generic Bessel component registration.
        ! params is evaluator-specific and interpreted by the selected
        ! Bessel application evaluator.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: component_index
        INTEGER, INTENT(IN) :: nparams
        REAL*8, INTENT(IN), DIMENSION(nparams) :: params

        IF (.NOT. COMPOSITE_BASIS_GRID_SET) STOP "initaxisymmetriccompositebasisexpansion must be called before addcompositebesselcomponent"
        IF (component_index < 1 .OR. component_index > COMPOSITE_NCOMP) STOP "invalid composite component_index"
        IF (nparams < 1) STOP "nparams must be >= 1 in addcompositebesselcomponent"
        IF (nparams > COMPOSITE_BESSEL_MAX_PARAMS) STOP "nparams exceeds COMPOSITE_BESSEL_MAX_PARAMS"

        COMPOSITE_BESSEL_PARAMS(:, component_index) = 0.0D0
        COMPOSITE_BESSEL_PARAMS(1:nparams, component_index) = params(1:nparams)
        COMPOSITE_BESSEL_NPARAMS(component_index) = nparams
        COMPOSITE_KIND(component_index) = BFE_KIND_BESSEL_DISK
        COMPOSITE_READY(component_index) = .TRUE.
        COMPOSITE_BASIS_FINALIZED = .FALSE.
    END SUBROUTINE addcompositebesselcomponent

    SUBROUTINE addcompositebesselexponentialdisk(component_index, sigma0, hR, hZ)
        ! Backward-compatible helper for the exponential-disk Bessel model.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: component_index
        REAL*8, INTENT(IN) :: sigma0, hR, hZ
        REAL*8, DIMENSION(4) :: params_disk

        IF (hR <= 0.0D0 .OR. hZ <= 0.0D0) STOP "hR and hZ must be positive"

        params_disk(1) = COMPOSITE_G
        params_disk(2) = sigma0
        params_disk(3) = hR
        params_disk(4) = hZ
        CALL addcompositebesselcomponent(component_index, params_disk, 4)
    END SUBROUTINE addcompositebesselexponentialdisk

    SUBROUTINE finalizeaxisymmetriccompositebasisexpansion()
        IMPLICIT NONE
        INTEGER :: i

        IF (.NOT. COMPOSITE_BASIS_GRID_SET) STOP "initaxisymmetriccompositebasisexpansion must be called before finalizeaxisymmetriccompositebasisexpansion"
        DO i = 1, COMPOSITE_NCOMP
            IF (.NOT. COMPOSITE_READY(i)) STOP "all composite components must be configured before finalizeaxisymmetriccompositebasisexpansion"
        END DO
        COMPOSITE_BASIS_FINALIZED = .TRUE.
    END SUBROUTINE finalizeaxisymmetriccompositebasisexpansion

    SUBROUTINE finalizecompositebasisexpansion()
        ! Backward-compatible alias.
        IMPLICIT NONE
        CALL finalizeaxisymmetriccompositebasisexpansion()
    END SUBROUTINE finalizecompositebasisexpansion

    SUBROUTINE axisymmetriccompositebasispotential_dispatch(params, N, x, y, z, ax, ay, az, phi)
        ! Dispatch wrapper that preserves the standard potential signature
        ! used by integrator procedure-pointer wiring.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(*) :: params
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi

        ! Keep params consumed so compilers do not warn in this wrapper.
        IF (params(1) /= params(1)) STOP "invalid NaN parameter in axisymmetriccompositebasispotential_dispatch"
        CALL axisymmetriccompositebasispotential(N, x, y, z, ax, ay, az, phi)
    END SUBROUTINE axisymmetriccompositebasispotential_dispatch

    SUBROUTINE compositebasispotential_dispatch(params, N, x, y, z, ax, ay, az, phi)
        ! Backward-compatible alias.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(*) :: params
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi

        CALL axisymmetriccompositebasispotential_dispatch(params, N, x, y, z, ax, ay, az, phi)
    END SUBROUTINE compositebasispotential_dispatch

    SUBROUTINE axisymmetriccompositebasispotential(N, x, y, z, ax, ay, az, phi)
        ! Evaluate all configured composite basis components and sum their forces.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi
        REAL*8, ALLOCATABLE :: ax_comp(:), ay_comp(:), az_comp(:), phi_comp(:)
        INTEGER :: i, nparams_bessel

        ax = 0.0D0
        ay = 0.0D0
        az = 0.0D0
        phi = 0.0D0

        IF (.NOT. COMPOSITE_BASIS_FINALIZED) STOP "finalizeaxisymmetriccompositebasisexpansion must be called before axisymmetriccompositebasispotential"

        ALLOCATE(ax_comp(N), ay_comp(N), az_comp(N), phi_comp(N))
        DO i = 1, COMPOSITE_NCOMP
            ax_comp = 0.0D0
            ay_comp = 0.0D0
            az_comp = 0.0D0
            phi_comp = 0.0D0

            SELECT CASE (COMPOSITE_KIND(i))
            CASE (BFE_KIND_LEGENDRE)
                CALL axisymmetricbasisexpansion_eval_component(N, x, y, z, COMPOSITE_R_GRID, &
                    COMPOSITE_PHI_L_GRID(:,:,i), COMPOSITE_DPHI_L_DR_GRID(:,:,i), &
                    ax_comp, ay_comp, az_comp, phi_comp)
            CASE (BFE_KIND_BESSEL_DISK)
                nparams_bessel = COMPOSITE_BESSEL_NPARAMS(i)
                IF (nparams_bessel < 1) STOP "invalid bessel parameter count in axisymmetriccompositebasispotential"
                CALL bessel_eval_component(exponential_disk_bessel_eval_component, &
                    COMPOSITE_BESSEL_PARAMS(1:nparams_bessel, i), N, &
                    x, y, z, ax_comp, ay_comp, az_comp, phi_comp)
            CASE DEFAULT
                STOP "unknown component kind in axisymmetriccompositebasispotential"
            END SELECT

            ax = ax + ax_comp
            ay = ay + ay_comp
            az = az + az_comp
            phi = phi + phi_comp
        END DO
        DEALLOCATE(ax_comp, ay_comp, az_comp, phi_comp)
    END SUBROUTINE axisymmetriccompositebasispotential

    SUBROUTINE compositebasispotential(N, x, y, z, ax, ay, az, phi)
        ! Backward-compatible alias.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi

        CALL axisymmetriccompositebasispotential(N, x, y, z, ax, ay, az, phi)
    END SUBROUTINE compositebasispotential

    SUBROUTINE bessel_eval_component(component_evaluator, params, N, x, y, z, ax, ay, az, phi)
        ! Generic Bessel-component dispatcher.
        ! Method: Bessel/Hankel representation used to solve Poisson for
        ! flattened axisymmetric components.
        ! Application: the specific density profile is implemented by the
        ! passed component_evaluator.
        IMPLICIT NONE
        PROCEDURE(bessel_component_evaluator_interface) :: component_evaluator
        REAL*8, INTENT(IN), DIMENSION(*) :: params
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi

        CALL component_evaluator(params, N, x, y, z, ax, ay, az, phi)
    END SUBROUTINE bessel_eval_component

    SUBROUTINE hernquist(params,N,x,y,z,ax,ay,az,phi)
        ! Hernquist potential
        ! params = [G, M, a]
        ! G = gravitational constant
        ! M = mass
        ! a = scale length
        ! x,y,z = coordinates
        ! ax,ay,az = acceleration

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8,INTENT(IN), DIMENSION(N) :: x,y,z
        REAL*8,INTENT(IN),dimension(3) :: params
        REAL*8,INTENT(OUT),DIMENSION(N) :: ax,ay,az,phi
        REAL*8, DIMENSION(N) :: r,amod
        REAL*8 :: G,M,a
        G = Params(1)
        M = params(2)
        a = params(3)

        r = sqrt(x*x + y*y + z*z)
        amod = -G*M / (r + a)**2


        ax = amod*x
        ay = amod*y
        az = amod*z
        phi = -G*M / (r + a)
    END SUBROUTINE hernquist
    
    

    SUBROUTINE plummer(params,N,x,y,z,ax,ay,az,phi)
        ! Plummer potential
        ! params = [G, M, a]
        ! G = gravitational constant
        ! M = mass
        ! b = scale length
        ! x,y,z = coordinates
        ! ax,ay,az = acceleration
        ! N = number of particles
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8,INTENT(IN), DIMENSION(N) :: x,y,z
        REAL*8,INTENT(IN),dimension(3) :: params
        REAL*8,INTENT(OUT),DIMENSION(N) :: ax,ay,az,phi
        REAL*8, DIMENSION(N) :: r,amod
        REAL*8 :: G,M,b
        G = Params(1)
        M = params(2) 
        b = params(3)

        r = sqrt(x*x + y*y + z*z)
        amod = -G*M / (r*r + b*b)**1.5
        
        ax = amod*x
        ay = amod*y
        az = amod*z
        phi = -G*M / (r*r + b*b)**0.5
    END SUBROUTINE plummer
    
    SUBROUTINE longmuralibar(params,N,x,y,z,ax,ay,az,phi)
        IMPLICIT NONE 
        INTEGER, INTENT(IN) :: N
        REAL*8,INTENT(IN), DIMENSION(N) :: x,y,z
        REAL*8,INTENT(IN),DIMENSION(5) :: params
        REAL*8,INTENT(OUT),DIMENSION(N) :: ax,ay,az,phi
        
        REAL*8 :: G,M,abar,bbar,cbar
        REAL*8,DIMENSION(N) :: Tplus,Tminus
        
        G=params(1)
        M=params(2)
        abar=params(3)
        bbar=params(4)
        cbar=params(5)

        Tplus=sqrt((abar+x)**2.+y*y+(bbar+sqrt(cbar*cbar+z*z))**2.)
        Tminus=sqrt((abar-x)**2.+y*y+(bbar+sqrt(cbar*cbar+z*z))**2.)
        phi=(G*M/2./abar)*log((x-abar+Tminus)/(x+abar+Tplus))


        ax=-2.*G*M*x/((Tplus*Tminus)*(Tplus+Tminus))
        ay=-G*M*y/((2.*Tplus*Tminus)*(y*y+(bbar+sqrt(z*z+cbar*cbar))**2.))*(Tplus+Tminus-4*x*x/(Tplus+Tminus))
        az=-G*M*z/((2.*Tplus*Tminus)*(y*y+(bbar+sqrt(z*z+cbar*cbar))**2.))*(Tplus+Tminus-4*x*x/(Tplus+Tminus))*&
        ((bbar+sqrt(z*z+cbar*cbar))/sqrt(z*z+cbar*cbar))


    END SUBROUTINE longmuralibar


    SUBROUTINE allensantillianhalo(params,N,x,y,z,ax,ay,az,phi)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8,INTENT(IN), DIMENSION(N) :: x,y,z
        REAL*8,INTENT(IN),dimension(5) :: params
        REAL*8,INTENT(OUT),DIMENSION(N) :: ax,ay,az,phi
        REAL*8, DIMENSION(N) :: term1
        REAL*8:: G,M,scale_length,exp,cutoffradius, Mtot
        REAL*8, DIMENSION(N) :: r,amod,d
        REAL*8:: term2,dcut
        LOGICAL, DIMENSION(N) :: outside_cutoff, at_zero

        REAL*8, DIMENSION(N) :: d_exp_minus_1, d_exp_minus_3

        G = params(1) ! gravitational constant
        M = params(2) ! mass parameter NOT total mass 
        scale_length = params(3) ! size parameter
        exp = params(4) ! exponential profile  (intended to be: 2.02)
        cutoffradius = params(5) ! cutoff radius (intended to be: 100 kpc)

        ! make dimensionless distance
        r = sqrt(x*x + y*y + z*z)
        d = r/scale_length
        dcut = cutoffradius/scale_length
        ! compute once for a speed up 
        d_exp_minus_1 = d**(exp-1)
        d_exp_minus_3 = d**(exp-3)
        
        Mtot = M * dcut**exp / (1 + dcut**(exp-1)) 
        amod = -(G*M/scale_length**3) * (d_exp_minus_3 / (1 + d_exp_minus_1))
        
        term1 = 1 + d_exp_minus_1
        term2 = 1 + dcut**(exp-1)

        phi = G*m/(scale_length*(exp-1)) * &
            log(term1/term2)  &
            - G*Mtot/cutoffradius

        ! Create masks for special cases
        outside_cutoff = (r > cutoffradius)
        at_zero = (r == 0)
        ! Handle special cases using MERGE instead of WHERE
        ! For points outside cutoff radius
        amod = MERGE(-(G*Mtot/r**3), amod, outside_cutoff)
        phi = MERGE(-G*Mtot/r, phi, outside_cutoff)

        ! For points at r=0 (prevent division by zero)
        amod = MERGE(0.0d0, amod, at_zero)
    
        ax=amod*x
        ay=amod*y
        az=amod*z

    end SUBROUTINE allensantillianhalo


    SUBROUTINE testallen(G,axhalo,Md_halo,N,x,y,z,ax,ay,az,phi)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8,INTENT(IN), DIMENSION(N) :: x,y,z
        REAL*8,INTENT(OUT),DIMENSION(N) :: ax,ay,az,phi
        REAL*8, DIMENSION(N) :: term1,term2
        REAL*8:: G,Md_halo,axhalo
        REAL*8, DIMENSION(N) :: r,Mhalo_r,amod,r3
        
        r = sqrt(x*x + y*y + z*z)
        r3=r*r*r

        Mhalo_r = Md_halo*(r/axhalo)**(2.02) ! the mass interior to r
        Mhalo_r = Mhalo_r/(1 + (r/axhalo)**(1.02))
        amod = -G*Mhalo_r/r3
        ax=amod*x
        ay=amod*y
        az=amod*z

        term1=-1.02/(1+(100./axhalo)**(1.02)) + log(1+(100./axhalo)**(1.02))
        term2=-1.02/(1+(r/axhalo)**(1.02)) + log(1+(r/axhalo)**(1.02))
        ! the equation derived from pouliasis et al 2017 and wrong allen santilian is weird
        ! deriving the potential from allen and martos 1986 makes more sense
        ! not clear to me from the articles why they have such a weird halo
        phi = -(G*Mhalo_r/r) - G*Md_halo/1.02/axhalo*(term1-term2)
    END SUBROUTINE testallen

    SUBROUTINE miyamotonagai(params,N,x,y,z,ax,ay,az,phi)
        IMPLICIT NONE
        REAL*8,INTENT(IN), DIMENSION(4) :: params
        INTEGER, INTENT(IN) :: N
        REAL*8,INTENT(IN), DIMENSION(N) :: x,y,z
        REAL*8,INTENT(OUT),DIMENSION(N) :: ax,ay,az,phi
        REAL*8, DIMENSION(N) :: R,amod,zmod
        REAL*8 :: G,M,a,b
        G = params(1) ! gravitational constant
        M = params(2) ! disk mass
        a = params(3) ! disk length
        b = params(4) ! disk height (thus, for galaxies a>b)
    
        R = sqrt(x*x + y*y)
        zmod =  a+(z*z + b*b)**0.5
        amod = -G*M / (R*R + zmod*zmod)**1.5
        ax=amod*x
        ay=amod*y
        az=amod*z*zmod/(sqrt(z*z+b*b))
        phi = -G*M / (R*R + zmod*zmod)**0.5
    END SUBROUTINE miyamotonagai


    SUBROUTINE pouliasis2017pii(params,N,x,y,z,ax,ay,az,phi)
        ! from pouliasis et al 2017
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8,INTENT(IN), DIMENSION(N) :: x,y,z
        REAL*8,INTENT(OUT),DIMENSION(N) :: ax,ay,az,phi
        REAL*8, INTENT(IN),DIMENSION(11) :: params
        REAL*8, DIMENSION(4) :: thindisk,thickdisk
        REAL*8, DIMENSION(5) :: halo
        REAL*8, DIMENSION(N) :: axH,ayH,azH,axD1,ayD1,azD1,axD2,ayD2,azD2
        REAL*8, DIMENSION(N) :: phiD1,phiD2,phiH
        ! HERE IS THE ORDER OF THE PARAMETERS
        ! params = [1, 2halo, 3halo, 4xp, 5utoffradius, 6disk1, 7calelength, 8caleheight, 9disk2, 10alelength, 11aleheight]        
        ! params = [G, Mhalo, ahalo, exp, cutoffradius, Mdisk1, scalelength, scaleheight, Mdisk2, scalelength, scaleheight]
        halo = (/params(1),params(2),params(3),params(4),params(5)/)
        thindisk = (/params(1),params(6),params(7),params(8)/)
        thickdisk = (/params(1),params(9),params(10),params(11)/)
        CALL allensantillianhalo(halo,N,x,y,z,axH,ayH,azH,phiH)
        CALL miyamotonagai(thindisk,N,x,y,z,axD1,ayD1,azD1,phiD1)
        CALL miyamotonagai(thickdisk,N,x,y,z,axD2,ayD2,azD2,phiD2)
        ax=axH+axD1+axD2
        ay=ayH+ayD1+ayD2
        az=azH+azD1+azD2
        phi=phiH+phiD1+phiD2

    END SUBROUTINE pouliasis2017pii

    SUBROUTINE NBODYPLUMMERS(params,N,x,y,z,ax,ay,az,phiTensor)
        ! Computeres the inter gravitational forces between N particles
        ! uses a plummer sphere for each particle. i.e. individual softening parameters for each body... 
        IMPLICIT NONE 
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN),DIMENSION(2*N + 1) :: params  ! first is the gravitational constant, then the masses, then the radii
        REAL*8, INTENT(IN),DIMENSION(N) :: x,y,z
        REAL*8, INTENT(OUT),DIMENSION(N) :: ax,ay,az
        REAL*8, INTENT(OUT),DIMENSION(N,N) :: phiTensor
        REAL*8, DIMENSION(N,N) :: FX,FY,FZ ! the forces on each particle
        REAL*8, DIMENSION(N) :: masses,scaleradii 
        REAL*8, DIMENSION(3) :: params2
        integer :: i,j
        masses=params(1+1:N+1)
        scaleradii=params(N+1+1:2*N+1) ! plus ones for the gravitational constant 
        ! force from i on all the others
        DO j=1,N
            params2(1)=params(1)
            params2(2)=masses(j)
            params2(3)=scaleradii(j)
            ! potential calculates the force of i on all the other particles
            CALL Plummer(params2,N,x-x(j),y-y(j),z-z(j),ax,ay,az,phiTensor(:,j))
            FX(:,j)=ax*masses
            FY(:,j)=ay*masses
            FZ(:,j)=az*masses
            phiTensor(j,j)=0
            FX(j,j)=0
            FY(j,j)=0
            FZ(j,j)=0
        END DO
        ! now sum up the forces
        ! we don't add a minus sign because we sum along the rows, which means force ON i from the others... 
        DO i=1,N
            ax(i)=sum(FX(i,:))/masses(i)
            ay(i)=sum(FY(i,:))/masses(i)
            az(i)=sum(FZ(i,:))/masses(i)
        END DO
    END SUBROUTINE NBODYPLUMMERS   
    
    SUBROUTINE pointmassconfiguration(G,NParticles,masses,xGC,yGC,zGC,N,x,y,z,ax,ay,az,phi)
        ! given a configuration of point masses,
        ! find the acceleration and potential at a given point
        IMPLICIT NONE 
        REAL*8, INTENT(IN) :: G
        INTEGER, INTENT(IN) :: N,NParticles
        REAL*8, INTENT(IN),DIMENSION(N) :: x,y,z
        REAL*8, INTENT(OUT),DIMENSION(N) :: ax,ay,az
        REAL*8, INTENT(OUT),DIMENSION(N) :: phi
        REAL*8, INTENT(IN),DIMENSION(NParticles) :: xGC, yGC, zGC, masses
        REAL*8 :: dx,dy,dz,dr,dr3
        INTEGER :: i,j

        ! initialize at zero
        ax=0
        ay=0
        az=0
        DO i=1,N
            DO j=1,NParticles
                dx=xGC(j)-x(i)
                dy=yGC(j)-y(i)
                dz=zGC(j)-z(i)
                dr=sqrt(dx*dx+dy*dy+dz*dz)
                dr3=dr*dr*dr
                ax(i)=ax(i)+G*masses(j)*dx/dr3
                ay(i)=ay(i)+G*masses(j)*dy/dr3
                az(i)=az(i)+G*masses(j)*dz/dr3
                phi(i)=phi(i)-G*masses(j)/dr
            END DO
        END DO
    end SUBROUTINE pointmassconfiguration





    SUBROUTINE default_init_basis_expansion(G)
        ! Auto-initialize with sensible defaults when the user has not called
        ! initaxisymmetricbasisexpansion explicitly.
        ! Grid: 100 log-spaced points from 1e-4 to 1e3; lmax = 20.
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: G
        INTEGER, PARAMETER :: default_lmax = 20
        INTEGER, PARAMETER :: default_nr   = 100
        REAL*8,  PARAMETER :: default_rmin = 1.0D-4
        REAL*8,  PARAMETER :: default_rmax = 1.0D3
        REAL*8 :: r_grid(default_nr)
        REAL*8 :: log_rmin, log_rmax, dlog_r
        INTEGER :: i

        log_rmin = LOG(default_rmin)
        log_rmax = LOG(default_rmax)
        dlog_r   = (log_rmax - log_rmin) / DBLE(default_nr - 1)
        DO i = 1, default_nr
            r_grid(i) = EXP(log_rmin + (i-1) * dlog_r)
        END DO
        CALL initaxisymmetricbasisexpansion(G, default_lmax, default_nr, r_grid)
    END SUBROUTINE default_init_basis_expansion

    SUBROUTINE project_exponential_oblate_halo(rho0, s0, q)
        ! Project rho(r,mu) = rho0*exp(-r/s0*sqrt(1-(1-1/q^2)*mu^2)) onto
        ! Legendre modes using Gauss-Legendre quadrature in mu.
        ! Fills BASIS_RHO_L_GRID for even l only.
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: rho0, s0, q
        INTEGER :: n_mu, i_r, l, k
        REAL*8, ALLOCATABLE :: mu_q(:), w_q(:), p(:)
        REAL*8 :: mu, rho_val, eta, factor
        REAL*8, PARAMETER :: pi_proj = 3.14159265358979323846D0

        n_mu = MAX(4 * (BASIS_LMAX + 1), 40)
        ALLOCATE(mu_q(n_mu), w_q(n_mu), p(0:BASIS_LMAX))
        CALL gauss_legendre_nodes_weights(n_mu, mu_q, w_q)

        ! eta = 1 - 1/q^2  so the density reads exp(-r/s0 * sqrt(1 - eta*mu^2))
        eta = 1.0D0 - 1.0D0 / (q * q)

        BASIS_RHO_L_GRID = 0.0D0
        DO i_r = 1, BASIS_NR
            DO k = 1, n_mu
                mu = mu_q(k)
                rho_val = rho0 * EXP(-BASIS_R_GRID(i_r) / s0 * SQRT(MAX(1.0D0 - eta*mu*mu, 0.0D0)))
                CALL legendre_p_all_axisymmetric(BASIS_LMAX, mu, p)
                DO l = 0, BASIS_LMAX, 2
                    factor = (2*l + 1) * 0.5D0 * w_q(k)
                    BASIS_RHO_L_GRID(l, i_r) = BASIS_RHO_L_GRID(l, i_r) + factor * rho_val * p(l)
                END DO
            END DO
        END DO

        DEALLOCATE(mu_q, w_q, p)
    END SUBROUTINE project_exponential_oblate_halo

    SUBROUTINE compute_phi_tables_from_rho()
        ! Build potential tables BASIS_PHI_L_GRID and BASIS_DPHI_L_DR_GRID from
        ! the projected density coefficients BASIS_RHO_L_GRID via the Green's
        ! function integrals for each l-mode:
        !   Phi_l(r) = -4*pi*G/(2l+1) * [ r^{-(l+1)} I_l^<(r) + r^l I_l^>(r) ]
        ! where I_l^<(r) = int_0^r  r'^{l+2}   rho_l(r') dr'
        !       I_l^>(r) = int_r^inf r'^{1-l} rho_l(r') dr'
        ! dPhi_l/dr = -4*pi*G/(2l+1) * [ -(l+1)*r^{-(l+2)} I_l^< + l*r^{l-1} I_l^> ]
        IMPLICIT NONE
        REAL*8, PARAMETER :: pi_phi = 3.14159265358979323846D0
        INTEGER :: l, i
        REAL*8  :: prefactor, r_lo, r_hi, f_lo, f_hi, dr
        REAL*8, ALLOCATABLE :: I_less(:), I_greater(:)

        ALLOCATE(I_less(BASIS_NR), I_greater(BASIS_NR))

        DO l = 0, BASIS_LMAX, 2
            prefactor = -4.0D0 * pi_phi * BASIS_G / DBLE(2*l + 1)

            ! Forward cumulative trapezoid: I_less(i) = int_0^r r'^{l+2} rho_l dr'
            I_less(1) = 0.0D0
            DO i = 2, BASIS_NR
                r_lo = BASIS_R_GRID(i-1);  r_hi = BASIS_R_GRID(i)
                f_lo = r_lo**(l+2) * BASIS_RHO_L_GRID(l, i-1)
                f_hi = r_hi**(l+2) * BASIS_RHO_L_GRID(l, i)
                I_less(i) = I_less(i-1) + 0.5D0 * (r_hi - r_lo) * (f_lo + f_hi)
            END DO

            ! Backward cumulative trapezoid: I_greater(i) = int_r^inf r'^{1-l} rho_l dr'
            I_greater(BASIS_NR) = 0.0D0
            DO i = BASIS_NR-1, 1, -1
                r_lo = BASIS_R_GRID(i);  r_hi = BASIS_R_GRID(i+1)
                f_lo = r_lo**(1-l) * BASIS_RHO_L_GRID(l, i)
                f_hi = r_hi**(1-l) * BASIS_RHO_L_GRID(l, i+1)
                dr = r_hi - r_lo
                I_greater(i) = I_greater(i+1) + 0.5D0 * dr * (f_lo + f_hi)
            END DO

            DO i = 1, BASIS_NR
                BASIS_PHI_L_GRID(l, i) = prefactor * ( &
                    BASIS_R_GRID(i)**(-(l+1)) * I_less(i) + &
                    BASIS_R_GRID(i)**l        * I_greater(i) )
                BASIS_DPHI_L_DR_GRID(l, i) = prefactor * ( &
                    -(l+1) * BASIS_R_GRID(i)**(-(l+2)) * I_less(i) + &
                    DBLE(l) * BASIS_R_GRID(i)**(l-1)   * I_greater(i) )
            END DO
        END DO

        DEALLOCATE(I_less, I_greater)
    END SUBROUTINE compute_phi_tables_from_rho

    SUBROUTINE axisymmetricbasisexpansion_eval(N, x, y, z, ax, ay, az, phi_out)
        ! Evaluate accelerations and potential for all N particles by
        ! interpolating the pre-computed l-mode tables and summing over modes.
        ! Cartesian force: a_i = -d_i phi, chain rule via (r, mu=z/r).
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN),  DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi_out
        INTEGER :: i, l, jlo, jhi, jmid
        REAL*8  :: r, r_safe, mu, R_cyl_sq
        REAL*8  :: alpha_r, phi_l_r, dphi_l_dr_r
        REAL*8  :: phi_v, dphi_dr_v, dphi_dmu_v
        REAL*8  :: p(0:BASIS_LMAX), dp_dmu(0:BASIS_LMAX)
        REAL*8, PARAMETER :: eps_r = 1.0D-30

        DO i = 1, N
            r      = SQRT(x(i)**2 + y(i)**2 + z(i)**2)
            r_safe = MAX(r, eps_r)
            mu     = z(i) / r_safe
            R_cyl_sq = x(i)**2 + y(i)**2

            ! Binary search for bracketing index in r_grid
            IF (r_safe <= BASIS_R_GRID(1)) THEN
                jlo = 1;  alpha_r = 0.0D0
            ELSE IF (r_safe >= BASIS_R_GRID(BASIS_NR)) THEN
                jlo = BASIS_NR - 1;  alpha_r = 1.0D0
            ELSE
                jlo = 1;  jhi = BASIS_NR
                DO WHILE (jhi - jlo > 1)
                    jmid = (jlo + jhi) / 2
                    IF (BASIS_R_GRID(jmid) <= r_safe) THEN
                        jlo = jmid
                    ELSE
                        jhi = jmid
                    END IF
                END DO
                ! Log-r interpolation weight (accurate on log-spaced grids)
                alpha_r = LOG(r_safe / BASIS_R_GRID(jlo)) / &
                          LOG(BASIS_R_GRID(jlo+1) / BASIS_R_GRID(jlo))
            END IF

            CALL legendre_axisymmetric_basis(BASIS_LMAX, mu, p, dp_dmu)

            phi_v = 0.0D0;  dphi_dr_v = 0.0D0;  dphi_dmu_v = 0.0D0
            DO l = 0, BASIS_LMAX, 2
                phi_l_r     = linear_interp_scalar(BASIS_PHI_L_GRID(l,jlo),     BASIS_PHI_L_GRID(l,jlo+1),     alpha_r)
                dphi_l_dr_r = linear_interp_scalar(BASIS_DPHI_L_DR_GRID(l,jlo), BASIS_DPHI_L_DR_GRID(l,jlo+1), alpha_r)
                phi_v      = phi_v      + phi_l_r     * p(l)
                dphi_dr_v  = dphi_dr_v  + dphi_l_dr_r * p(l)
                dphi_dmu_v = dphi_dmu_v + phi_l_r     * dp_dmu(l)
            END DO

            phi_out(i) = phi_v
            ! a = -grad(phi).  With phi(r,mu), mu=z/r:
            !   d_phi/dx = phi_r*(x/r) - phi_mu*(z*x/r^3)
            !   d_phi/dy = phi_r*(y/r) - phi_mu*(z*y/r^3)
            !   d_phi/dz = phi_r*(z/r) + phi_mu*(R^2/r^3)
            ax(i) = -dphi_dr_v * x(i)/r_safe + dphi_dmu_v * z(i)*x(i)/r_safe**3
            ay(i) = -dphi_dr_v * y(i)/r_safe + dphi_dmu_v * z(i)*y(i)/r_safe**3
            az(i) = -dphi_dr_v * z(i)/r_safe - dphi_dmu_v * R_cyl_sq/r_safe**3
        END DO
    END SUBROUTINE axisymmetricbasisexpansion_eval

    SUBROUTINE compute_phi_tables_from_rho_component(r_grid, G, rho_l_grid, phi_l_grid, dphi_l_dr_grid)
        IMPLICIT NONE
        REAL*8, INTENT(IN), DIMENSION(:) :: r_grid
        REAL*8, INTENT(IN) :: G
        REAL*8, INTENT(IN), DIMENSION(:,:) :: rho_l_grid
        REAL*8, INTENT(OUT), DIMENSION(:,:) :: phi_l_grid, dphi_l_dr_grid
        REAL*8, PARAMETER :: pi_phi = 3.14159265358979323846D0
        INTEGER :: l, i, l_lo, l_hi, nr
        REAL*8  :: prefactor, r_lo, r_hi, f_lo, f_hi, dr
        REAL*8, ALLOCATABLE :: I_less(:), I_greater(:)

        nr = SIZE(r_grid)
        l_lo = LBOUND(rho_l_grid, 1)
        l_hi = UBOUND(rho_l_grid, 1)

        ALLOCATE(I_less(nr), I_greater(nr))

        DO l = l_lo, l_hi, 2
            prefactor = -4.0D0 * pi_phi * G / DBLE(2*l + 1)

            I_less(1) = 0.0D0
            DO i = 2, nr
                r_lo = r_grid(i-1);  r_hi = r_grid(i)
                f_lo = r_lo**(l+2) * rho_l_grid(l, i-1)
                f_hi = r_hi**(l+2) * rho_l_grid(l, i)
                I_less(i) = I_less(i-1) + 0.5D0 * (r_hi - r_lo) * (f_lo + f_hi)
            END DO

            I_greater(nr) = 0.0D0
            DO i = nr-1, 1, -1
                r_lo = r_grid(i);  r_hi = r_grid(i+1)
                f_lo = r_lo**(1-l) * rho_l_grid(l, i)
                f_hi = r_hi**(1-l) * rho_l_grid(l, i+1)
                dr = r_hi - r_lo
                I_greater(i) = I_greater(i+1) + 0.5D0 * dr * (f_lo + f_hi)
            END DO

            DO i = 1, nr
                phi_l_grid(l, i) = prefactor * ( &
                    r_grid(i)**(-(l+1)) * I_less(i) + &
                    r_grid(i)**l        * I_greater(i) )
                dphi_l_dr_grid(l, i) = prefactor * ( &
                    -(l+1) * r_grid(i)**(-(l+2)) * I_less(i) + &
                    DBLE(l) * r_grid(i)**(l-1)   * I_greater(i) )
            END DO
        END DO

        DEALLOCATE(I_less, I_greater)
    END SUBROUTINE compute_phi_tables_from_rho_component

    SUBROUTINE axisymmetricbasisexpansion_eval_component(N, x, y, z, r_grid, phi_l_grid, dphi_l_dr_grid, ax, ay, az, phi_out)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN),  DIMENSION(N) :: x, y, z
        REAL*8, INTENT(IN),  DIMENSION(:) :: r_grid
        REAL*8, INTENT(IN),  DIMENSION(:,:) :: phi_l_grid, dphi_l_dr_grid
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi_out
        INTEGER :: i, l, jlo, jhi, jmid, l_lo, l_hi, nr
        REAL*8  :: r, r_safe, mu, R_cyl_sq
        REAL*8  :: alpha_r, phi_l_r, dphi_l_dr_r
        REAL*8  :: phi_v, dphi_dr_v, dphi_dmu_v
        REAL*8, ALLOCATABLE :: p(:), dp_dmu(:)
        REAL*8, PARAMETER :: eps_r = 1.0D-30

        nr = SIZE(r_grid)
        l_lo = LBOUND(phi_l_grid, 1)
        l_hi = UBOUND(phi_l_grid, 1)
        ALLOCATE(p(l_lo:l_hi), dp_dmu(l_lo:l_hi))

        DO i = 1, N
            r      = SQRT(x(i)**2 + y(i)**2 + z(i)**2)
            r_safe = MAX(r, eps_r)
            mu     = z(i) / r_safe
            R_cyl_sq = x(i)**2 + y(i)**2

            IF (r_safe <= r_grid(1)) THEN
                jlo = 1;  alpha_r = 0.0D0
            ELSE IF (r_safe >= r_grid(nr)) THEN
                jlo = nr - 1;  alpha_r = 1.0D0
            ELSE
                jlo = 1;  jhi = nr
                DO WHILE (jhi - jlo > 1)
                    jmid = (jlo + jhi) / 2
                    IF (r_grid(jmid) <= r_safe) THEN
                        jlo = jmid
                    ELSE
                        jhi = jmid
                    END IF
                END DO
                alpha_r = LOG(r_safe / r_grid(jlo)) / LOG(r_grid(jlo+1) / r_grid(jlo))
            END IF

            CALL legendre_axisymmetric_basis(l_hi, mu, p, dp_dmu)

            phi_v = 0.0D0;  dphi_dr_v = 0.0D0;  dphi_dmu_v = 0.0D0
            DO l = l_lo, l_hi, 2
                phi_l_r     = linear_interp_scalar(phi_l_grid(l,jlo),     phi_l_grid(l,jlo+1),     alpha_r)
                dphi_l_dr_r = linear_interp_scalar(dphi_l_dr_grid(l,jlo), dphi_l_dr_grid(l,jlo+1), alpha_r)
                phi_v      = phi_v      + phi_l_r     * p(l)
                dphi_dr_v  = dphi_dr_v  + dphi_l_dr_r * p(l)
                dphi_dmu_v = dphi_dmu_v + phi_l_r     * dp_dmu(l)
            END DO

            phi_out(i) = phi_v
            ax(i) = -dphi_dr_v * x(i)/r_safe + dphi_dmu_v * z(i)*x(i)/r_safe**3
            ay(i) = -dphi_dr_v * y(i)/r_safe + dphi_dmu_v * z(i)*y(i)/r_safe**3
            az(i) = -dphi_dr_v * z(i)/r_safe - dphi_dmu_v * R_cyl_sq/r_safe**3
        END DO

        DEALLOCATE(p, dp_dmu)
    END SUBROUTINE axisymmetricbasisexpansion_eval_component

    SUBROUTINE exponential_disk_bessel_eval_component(params, N, x, y, z, ax, ay, az, phi)
        ! Exponential-disk application of the generic Bessel/Poisson method.
        ! This is the preferred basis path for strongly flattened disk-like
        ! components (e.g. q <= 0.3), while Legendre is better for near-
        ! spherical components.
        !
        ! Density model associated with this evaluator:
        !   rho(R,z) = sigma0/(2*hZ) * exp(-R/hR - |z|/hZ)
        IMPLICIT NONE
        REAL*8, INTENT(IN), DIMENSION(*) :: params
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN), DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi
        INTEGER, PARAMETER :: nk = 256
        REAL*8, PARAMETER :: pi = 3.14159265358979323846D0
        REAL*8, PARAMETER :: eps = 1.0D-16
        INTEGER :: i, j
        REAL*8 :: R, absz, signz, aR, k, dk, kmax, A
        REAL*8 :: G, sigma0, hR, hZ
        REAL*8 :: sum_phi, sum_ar, sum_az, weight
        REAL*8 :: kernel, j0, j1
        REAL*8 :: scale_min

        G = params(1)
        sigma0 = params(2)
        hR = params(3)
        hZ = params(4)

        scale_min = MAX(MIN(hR, hZ), 1.0D-6)
        kmax = 60.0D0 / scale_min
        dk = kmax / DBLE(nk - 1)
        A = 2.0D0 * pi * G * sigma0 * hR * hR

        DO i = 1, N
            R = SQRT(x(i)**2 + y(i)**2)
            absz = ABS(z(i))
            IF (z(i) > 0.0D0) THEN
                signz = 1.0D0
            ELSE IF (z(i) < 0.0D0) THEN
                signz = -1.0D0
            ELSE
                signz = 0.0D0
            END IF

            sum_phi = 0.0D0
            sum_ar = 0.0D0
            sum_az = 0.0D0
            DO j = 1, nk
                k = DBLE(j-1) * dk
                IF (j == 1 .OR. j == nk) THEN
                    weight = 0.5D0
                ELSE
                    weight = 1.0D0
                END IF

                kernel = EXP(-k * absz) / ( (1.0D0 + (k*hR)**2)**1.5D0 * (1.0D0 + k*hZ) )
                j0 = bessel_j0_scalar(k * R)
                j1 = bessel_j1_scalar(k * R)

                sum_phi = sum_phi + weight * j0 * kernel
                sum_ar  = sum_ar  + weight * k * j1 * kernel
                sum_az  = sum_az  + weight * k * j0 * kernel
            END DO

            phi(i) = -A * dk * sum_phi
            aR = -A * dk * sum_ar
            az(i) = -A * dk * signz * sum_az

            IF (R > eps) THEN
                ax(i) = aR * x(i) / R
                ay(i) = aR * y(i) / R
            ELSE
                ax(i) = 0.0D0
                ay(i) = 0.0D0
            END IF
        END DO
    END SUBROUTINE exponential_disk_bessel_eval_component

    SUBROUTINE exponential_oblate_halo(params, N, x, y, z, ax, ay, az, phi)
        ! Axisymmetric exponential oblate halo:
        !   rho(R,z) = rho0 * exp(-1/s0 * sqrt(R^2 + z^2/q^2))
        ! params = [G, rho0, s0, q]
        !
        ! On first call (or after clear): if no grid has been set up, auto-
        ! initializes with defaults (lmax=20, 100 log-spaced points 1e-4..1e3).
        ! Then projects rho onto Legendre modes and computes potential tables.
        ! Subsequent calls skip projection and go straight to evaluation.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN),  DIMENSION(4) :: params
        REAL*8, INTENT(IN),  DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi
        REAL*8 :: G, rho0, s0, q

        G    = params(1)
        rho0 = params(2)
        s0   = params(3)
        q    = params(4)

        IF (.NOT. BASIS_EXPANSION_INITIALIZED) THEN
            IF (.NOT. BASIS_GRID_SET) THEN
                CALL default_init_basis_expansion(G)
            END IF
            CALL project_exponential_oblate_halo(rho0, s0, q)
            CALL compute_phi_tables_from_rho()
            BASIS_EXPANSION_INITIALIZED = .TRUE.
        END IF

        CALL axisymmetricbasisexpansion_eval(N, x, y, z, ax, ay, az, phi)
    END SUBROUTINE exponential_oblate_halo


    SUBROUTINE project_ibata2024halo(rho0, r0, rt, q, gamma, beta)
        ! Project
        ! rho(s) = rho0 * (s/r0)^(-gamma) * (1 + s/r0)^(gamma-beta) * exp(-(s/rt)^2)
        ! with s = r * sqrt(1 - eta*mu^2), eta = 1 - 1/q^2
        ! onto even-l Legendre modes, filling BASIS_RHO_L_GRID.
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: rho0, r0, rt, q, gamma, beta
        INTEGER :: n_mu, i_r, l, k
        REAL*8, ALLOCATABLE :: mu_q(:), w_q(:), p(:)
        REAL*8 :: mu, rho_val, eta, factor, s, x
        REAL*8, PARAMETER :: s_floor = 1.0D-12

        n_mu = MAX(4 * (BASIS_LMAX + 1), 40)
        ALLOCATE(mu_q(n_mu), w_q(n_mu), p(0:BASIS_LMAX))
        CALL gauss_legendre_nodes_weights(n_mu, mu_q, w_q)

        eta = 1.0D0 - 1.0D0 / (q * q)

        BASIS_RHO_L_GRID = 0.0D0
        DO i_r = 1, BASIS_NR
            DO k = 1, n_mu
                mu = mu_q(k)

                s = BASIS_R_GRID(i_r) * SQRT(MAX(1.0D0 - eta*mu*mu, 0.0D0))
                x = MAX(s / r0, s_floor)

                rho_val = rho0 * x**(-gamma) * (1.0D0 + x)**(gamma - beta) * EXP(-(s/rt)**2)

                CALL legendre_p_all_axisymmetric(BASIS_LMAX, mu, p)
                DO l = 0, BASIS_LMAX, 2
                    factor = (2*l + 1) * 0.5D0 * w_q(k)
                    BASIS_RHO_L_GRID(l, i_r) = BASIS_RHO_L_GRID(l, i_r) + factor * rho_val * p(l)
                END DO
            END DO
        END DO

        DEALLOCATE(mu_q, w_q, p)
    END SUBROUTINE project_ibata2024halo    

    SUBROUTINE ibata2024halo(params, N, x, y, z, ax, ay, az, phi)
        ! Axisymmetric Ibata-like halo:
        ! rho(s) = rho0 * (s/r0)^(-gamma) * (1 + s/r0)^(gamma-beta) * exp(-(s/rt)^2)
        ! s = sqrt(R^2 + z^2/q^2), R^2 = x^2 + y^2
        !
        ! params = [G, rho0, r0, rt, q, gamma, beta]
        !
        ! Same initialization logic as exponential_oblate_halo.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN),  DIMENSION(7) :: params
        REAL*8, INTENT(IN),  DIMENSION(N) :: x, y, z
        REAL*8, INTENT(OUT), DIMENSION(N) :: ax, ay, az, phi
        REAL*8 :: G, rho0, r0, rt, q, gamma, beta

        G     = params(1)
        rho0  = params(2)
        r0    = params(3)
        rt    = params(4)
        q     = params(5)
        gamma = params(6)
        beta  = params(7)

        IF (.NOT. BASIS_EXPANSION_INITIALIZED) THEN
            IF (.NOT. BASIS_GRID_SET) THEN
                CALL default_init_basis_expansion(G)
            END IF
            CALL project_ibata2024halo(rho0, r0, rt, q, gamma, beta)
            CALL compute_phi_tables_from_rho()
            BASIS_EXPANSION_INITIALIZED = .TRUE.
        END IF

        CALL axisymmetricbasisexpansion_eval(N, x, y, z, ax, ay, az, phi)
    END SUBROUTINE ibata2024halo

end module potentials




