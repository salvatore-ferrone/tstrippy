import numpy as np


#############################################################
#################### GENERATE CONDITIONS ####################
#############################################################

def isotropicplummer(G,Mt,halfmassradius,NP):
    '''
    purpose: generate initial conditions for a Plummer sphere
    input:  Mt = total mass of the sphere
            halfmassradius = observed half mass radius (not plummer radius)
            NP = number of particles
            G = gravitational constant
    output: x,y,z,vx,vy,vz = initial conditions
    '''
    # random uniform sampling of the enclused mass
    rc=convertHalfMassRadiusToPlummerRadius(halfmassradius)
    Mxs=np.random.rand(NP)
    # inverse transform sampling for the radius
    r=PlummerRadius(Mxs,rc)
    phi,theta,_ = UniformSphere(NP)
    # generate x,y,z
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    vx,vy,vz,_=velocitySampling(DFPlummer,Mt,rc,r,G,storeCDF=False)
    return x,y,z,vx,vy,vz

#############################################################
##################### GENERAL FUNCTIONS #####################
#############################################################
def UniformSphere(NP):
    """
    Generate a uniform distribution of points on a sphere
    Parameter:
    NP: number of particles
    RETURNS
        phi: angle in longitude
        theta: angle in latitude
        thetaWRONG: incorrect latitude sampling
    """
    # uniform angle in longitude
    phi = np.random.rand(NP)*2*np.pi
    # inverse transform sampling for the angle in latitude
    theta = np.arccos(2*np.random.rand(NP)-1)
    # generate the incorrect latitude sampling for comparison
    thetaWRONG = np.pi*np.random.rand(NP) 
    return phi,theta,thetaWRONG

def EscapeVelocity(PotentialEnergy,Mt,a,r,G):
    """
    the escape velocity
    """
    # assert that PotentialEnergy is a function
    assert callable(PotentialEnergy), "PotentialEnergy is not a function"

    return np.sqrt(-2*PotentialEnergy(Mt,a,r,G))

#############################################################
##################### PLUMMER FUNCTIONS #####################
#############################################################
def DFPlummer(Mt,a,r,v,G):
    """
    The DF
    """
    E=(1/2)*v**2 + PlummerPotentialEnergy(Mt,a,r,G)
    return (-E)**(7/5)

def PlummerPotentialEnergy(Mt,a,r,G):
    """
    Mt: total mass
    a = plummer scale length
    r: radius
    G = gravitational constant

    """
    return -G*Mt/np.sqrt(r**2 + a**2)

def PlummerMassProfile(Mt,Rc,r):
    """
    the amount of enclosed mass as a function of r 
    Parameter:
    M: total mass of the Plummer sphere
    Rc: core radius of the Plummer sphere
    r: radius at which we want to know the enclosed mass
    """
    return Mt * r**3 * (r**2 + Rc**2 )**(-3/2)   

def PlummerRadius(Mx,Rc):
    """
    The radius at which we enclose Menc
    Parameter:
        Mx: the fraction of the total mass we want to enclose
        Rc: core radius of the Plummer sphere
    """
    return Rc*np.sqrt((Mx**(2/3))/(1-Mx**(2/3)))

def convertHalfMassRadiusToPlummerRadius(halfmassradius):
    """
    purpose: convert the half mass radius to the plummer radius
    input: halfmassradius = observed half mass radius (not plummer radius)
    output: a = plummer radius
    """
    const = ((1/(1/2)**(2/3))-1)**(1/2)
    a=const*halfmassradius
    return a

###############################################
######### VELOCITY SAMPLING FUNCTIONS #########
###############################################
def velocityCDF(vesc,DistFunc,Mt,Rc,r,G,relNumGaurd=1000,nsamp=1000):
    """
    purpose: generate the CDF for the velocity distribution function
    INPUTS:
        vesc: escape velocity at the radius r
        DistFunc: the distribution function
        Mt: total mass of the system
        Rc: core radius of the Plummer sphere
        r: radius of the particle of interest
        G: gravitational constant
        relNumGaurd: The numerical gaurd for not including vesc in the grid
    OUTPUTS:
        CDFnorm: the normalized CDF
        vels: the velocity grid
    """
    # assert that DF is a function
    assert callable(DistFunc), "DistFunc is not a function"
    dv = vesc/relNumGaurd
    vels=np.linspace(0,vesc-dv,nsamp)
    DF=DistFunc(Mt,Rc,r,vels,G)
    # remember, the velocity sphere, for equal probability
    CDF=np.cumsum(DF*(vels**2)*np.pi*4)
    CDFmax=CDF.max()
    CDFmin=CDF.min()
    CDFnorm=(CDF-CDFmin)/(CDFmax-CDFmin)
    return CDFnorm,vels
    
def velocitySampling(DistFunc,Mt,Rc,rAll,G,storeCDF=True):
    """ 
    Sample the velocities of a particle distribution
    INPUTS:
        DisFunc: the distribution function
        Mt: total mass of the system
        Rc: core radius of the Plummer sphere
        rAll: array of radii of the particles of interest
        G: gravitational constant
        storeCDF: boolean to store the CDF
    OUTPUTS:    
        vx,vy,vz: arrays of the sampled velocities
        OPTIONAL
            CDFAll: all cummulative distribution functions for each particle in rAll
            vels: all of the velocity grids for each particle in rAll
            rannums (np.array, same size as rAll): each random number corresponding to the random inverse transform sampling each CDF
    """
    
    if DistFunc.__name__=='DFPlummer':
        PotentialEnergy=PlummerPotentialEnergy
    else:
        raise ValueError('DistFunc not recognized')
    # number of samples to take for the CDF
    nsamp=1000
    # numerical gaurd for not including vesc in the grid
    relNumGaurd=1000    
    # get the number of particles
    NP=len(rAll)
    # get the v_escape for each particle
    vesc=EscapeVelocity(PotentialEnergy,Mt,Rc,rAll,G)
    # initialize velocity array
    speed=np.zeros(NP)

    if storeCDF:
        testVels=np.zeros((NP,nsamp))
        CDFAll=np.zeros((NP,nsamp))
        # store the random numbers
        rannums=np.zeros(NP)
        # loop over the particles
        for i in range(NP):
            CDFnorm,vels=velocityCDF(\
                vesc[i],DistFunc,Mt,Rc,rAll[i],G=G,relNumGaurd=relNumGaurd,nsamp=nsamp)
            # get the random number
            rannum=np.random.rand(1)
            # inverse transform sample the velocity
            speed[i]=np.interp(rannum,CDFnorm,vels)
            # store the CDF and velocity grid
            CDFAll[i,:]=CDFnorm
            testVels[i,:]=vels
            # store the random number
            rannums[i]=rannum
    else:
        # don't store the extra information
        for i in range(NP):
            CDFnorm,vels=velocityCDF(vesc[i],DistFunc,Mt,Rc,rAll[i],G=G,relNumGaurd=relNumGaurd,nsamp=nsamp)
            # get the random number
            rannum=np.random.rand(1)
            # inverse transform sample the velocity
            speed[i]=np.interp(rannum[0],CDFnorm,vels)
        
    # get the angles for the velocity
    vel=np.random.random((NP,3))-0.5
    vel=vel/np.sqrt(np.sum(vel**2,axis=1))[:,None]
    vel=vel*speed[:,None]
    vx,vy,vz=vel[:,0],vel[:,1],vel[:,2]

    if storeCDF:
        return vx,vy,vz,speed,CDFAll,testVels,rannums
    else:
        return vx,vy,vz,speed
