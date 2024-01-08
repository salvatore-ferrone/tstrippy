import astropy.units as u
from astropy.io import fits
import numpy as np
import pkg_resources


class baumgardtMWGCs:
    """
    A class to handle data from the Baumgardt Milky Way Globular Clusters (MWGCs).

    This class provides methods to extract data from a FITS file, select a globular cluster of interest, 
    and build a covariance matrix of the observed kinematics of a given globular cluster.

    Attributes
    ----------
    pathtoclusterdata : str
        The path to the FITS file containing the cluster data.
    data : dict
        The data extracted from the FITS file.
    units : dict
        The units for the data fields.
    GCdex : int
        The index of the selected globular cluster in the data.
    GCname : str
        The name of the selected globular cluster.
    RA : astropy.units.quantity.Quantity
        The right ascension of the selected globular cluster.
    ERRA : astropy.units.quantity.Quantity
        The error on the right ascension of the selected globular cluster.
    DEC : astropy.units.quantity.Quantity
        The declination of the selected globular cluster.
    ERDEC : astropy.units.quantity.Quantity
        The error on the declination of the selected globular cluster.
    Rsun : astropy.units.quantity.Quantity
        The distance from the Sun to the selected globular cluster.
    ERsun : astropy.units.quantity.Quantity
        The error on the distance from the Sun to the selected globular cluster.
    RV : astropy.units.quantity.Quantity
        The radial velocity of the selected globular cluster.
    ERV : astropy.units.quantity.Quantity
        The error on the radial velocity of the selected globular cluster.
    mualpha : astropy.units.quantity.Quantity
        The proper motion in right ascension of the selected globular cluster.
    ERmualpha : astropy.units.quantity.Quantity
        The error on the proper motion in right ascension of the selected globular cluster.
    mu_delta : astropy.units.quantity.Quantity
        The proper motion in declination of the selected globular cluster.
    ERmu_delta : astropy.units.quantity.Quantity
        The error on the proper motion in declination of the selected globular cluster.
    rhopmrade : astropy.units.quantity.Quantity
        The correlation between the proper motions in right ascension and declination of the selected globular cluster.
    Mass : astropy.units.quantity.Quantity
        The mass of the selected globular cluster.
    DM : astropy.units.quantity.Quantity
        The distance modulus of the selected globular cluster.
    rh_m : astropy.units.quantity.Quantity
        The half-mass radius of the selected globular cluster.
    ERrh_m : astropy.units.quantity.Quantity
        The error on the half-mass radius of the selected globular cluster.
    """

    def __init__(
        self,
        units ={
            "RA":u.degree,
            "DEC":u.degree,
            "Rsun":u.kpc,
            "RV":u.km/u.s,
            "mualpha":u.mas/u.yr,
            "mu_delta":u.mas/u.yr,
            "Mass":u.Msun,
            "rh_m":u.kpc,
            "DM":u.Msun,
            "rhopmrade":""}
        ):
        """ extracts the data from the fits file and stores it in a dictionary

        """        
        
        self.pathtoclusterdata=pkg_resources.resource_filename('tstrippy', 'data/2023-03-28-merged.fits')
        
        self._extractdatafromfits()
        self.units=units


    def _extractdatafromfits(self):
        """ internal function to extract the data from the fits file and store it in a dictionary
        """        
        GCfits=fits.open(self.pathtoclusterdata)
        fields = [GCfits[1].header['TTYPE'+str(i)] for i in range(1,GCfits[1].header['TFIELDS']+1)]
        units = [GCfits[1].header['TUNIT'+str(i)] for i in range(1,GCfits[1].header['TFIELDS']+1)]
        data={}
        for i in range(len(fields)):
            if units[i]=="Name":
                data[fields[i]] = GCfits[1].data[fields[i]]
            else:
                if units[i]=="unitless":
                    units[i]=""
                data[fields[i]] = GCfits[1].data[fields[i]]*u.Unit(units[i])
        self.data=data
        GCfits.close()
    
    def SetGlobularClusterOfInterest(self,GCname):
        """Pick a globular cluster of interest and extract its data to tbe class directly

        Parameters
        ----------
        GCname : str
            name of the globular cluter of interest. Must be one of the names in the fits file
        
        Returns
        -------
        None
            None
        """        

        self.GCdex=np.where(self.data['Cluster']==GCname)[0][0]
        self.GCname=GCname
        self.RA=self.data['RA'][self.GCdex].to(self.units['RA'])
        self.ERRA=0*self.units['RA'] # no errors on these
        self.DEC=self.data['DEC'][self.GCdex].to(self.units['DEC'])
        self.ERDEC=0*self.units['DEC'] # no errors on these
        self.Rsun=self.data['Rsun'][self.GCdex].to(self.units['Rsun'])
        self.ERsun=self.data['ERsun'][self.GCdex].to(self.units['Rsun'])
        self.RV=self.data['RV'][self.GCdex].to(self.units['RV'])
        self.ERV=self.data['ERV'][self.GCdex].to(self.units['RV'])
        self.mualpha=self.data['mualpha'][self.GCdex].to(self.units['mualpha'])
        self.ERmualpha=self.data['ERmualpha'][self.GCdex].to(self.units['mualpha'])
        self.mu_delta=self.data['mu_delta'][self.GCdex].to(self.units['mu_delta'])
        self.ERmu_delta=self.data['ERmu_delta'][self.GCdex].to(self.units['mu_delta'])
        self.rhopmrade=self.data['rhopmrade'][self.GCdex].to(self.units['rhopmrade'])
        self.Mass=self.data['Mass'][self.GCdex].to(self.units['Mass'])
        self.DM=self.data['DM'][self.GCdex].to(self.units['DM'])
        self.rh_m=self.data['rh_m'][self.GCdex].to(self.units['rh_m'])
        self.ERrh_m=0*self.rh_m.unit # no errors on these
    
    def BuildCovarianceMatrix(self,ERRA=0,ERDEC=0,ERhalfmassradius=0):
        """build the covariance matrix of the observed kinematics of a given GC. Intended to be used for sampling

        Parameters
        ----------
        ERRA : int, optional
            error on the Right Ascension, by default 0
        ERDEC : int, optional
            error on the Declination, by default 0
        ERhalfmassradius : int, optional
            Error of the half mass radius, by default 0
        
        Examples
        --------
        >>> self.BuildCovarianceMatrix()
        >>> RA,DEC,D,RV,mualpha,mu_delta,Mass,halfmassradius = np.random.multivariate_normal(self.meanParameters,self.covarianceParameters)


        Returns
        -------
        None
            None
            
        Note
        ----
        the order of the data are:
            [RA,DEC,Rsun,RV,mualpha,mu_delta,rh_m]
        the data can be accessed with the following:
            self.sigmaParameters = sigmas
            self.meanParameters = means
            self.covarianceParameters =

        """        
        means = [self.RA.to(self.units['RA']).value,
            self.DEC.to(self.units['DEC']).value,
            self.Rsun.to(self.units['Rsun']).value,
            self.RV.to(self.units['RV']).value,
            self.mualpha.to(self.units['mualpha']).value,
            self.mu_delta.to(self.units['mu_delta']).value,
            self.Mass.to(self.units['Mass']).value,
            self.rh_m.to(self.units['rh_m']).value,]
        sigmas = [self.ERRA.to(self.units['RA']).value,
            self.ERDEC.to(self.units['DEC']).value,
            self.ERsun.to(self.units['Rsun']).value,
            self.ERV.to(self.units['RV']).value,
            self.ERmualpha.to(self.units['mualpha']).value,
            self.ERmu_delta.to(self.units['mu_delta']).value,
            self.DM.to(self.units['Mass']).value,
            self.ERrh_m.to(self.units['rh_m']).value,]
        Nparams = len(means)
        # put in the variances
        covarianceMatrix=np.zeros((Nparams,Nparams))
        for i in range(Nparams):
            covarianceMatrix[i,i]=sigmas[i]**2
        # convert correlation to covariance and put in the covariance
        covarianceMatrix[4,5]=covarianceMatrix[5,4]=self.rhopmrade*self.ERmualpha.value*self.ERmu_delta.value # only the covariance between the proper motions
        self.fields = {
            "RA":{"name":"RA","unit":self.units['RA']},
            "DEC":{"name":"DEC","unit":self.units['DEC']},
            "Rsun":{"name":"Rsun","unit":self.units['Rsun']},
            "RV":{"name":"RV","unit":self.units['RV']},
            "mualpha":{"name":r"$\mu_{\alpha}cos{\delta}$","unit":self.units['mualpha']},
            "mu_delta":{"name":r"$\mu_{\delta}$","unit":self.units['mu_delta']},
            "Mass":{"name":"Mass","unit":r"$M_{\odot}$"},
            "rh_m":{"name":r"$r_{\frac{1}{2}M}$","unit":self.units['rh_m']},
        }
        self.sigmaParameters = sigmas
        self.meanParameters = means
        self.covarianceParameters = covarianceMatrix


