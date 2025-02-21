MODULE potentials
    REAL*8, PARAMETER :: pi = 2.0d0 * acos(0.0d0)
    contains
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
        REAL*8, DIMENSION(N) :: term1,term2
        REAL*8:: G,M,a,exp,cutoffradius
        REAL*8, DIMENSION(N) :: r,massMod,amod
        G = params(1) ! gravitational constant
        M = params(2) ! total halo mass
        a = params(3) ! size parameter
        exp = params(4) ! exponential profile  (intended to be: 2.02)
        cutoffradius = params(5) ! cutoff radius (intended to be: 100 kpc)

        r = sqrt(x*x + y*y + z*z)
        massMod = M*(r/a)**(exp) ! the mass interior to r
        massMod = massMod/(1 + (r/a)**(exp-1))
        amod = -G*massMod/r**3
        ax=amod*x
        ay=amod*y
        az=amod*z

        term1=-(exp-1)/(1+(cutoffradius/a)**(exp-1)) + log(1+(cutoffradius/a)**(exp-1))
        term2=-(exp-1)/(1+(r/a)**(exp-1)) + log(1+(r/a)**(exp-1))
        ! the equation derived from pouliasis et al 2017 and wrong allen santilian is weird
        ! deriving the potential from allen and martos 1986 makes more sense
        ! not clear to me from the articles why they have such a weird halo
        phi = -(G*massMod/r) - (G*M/((exp-1)*a))*(term1-term2)
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





end module potentials




