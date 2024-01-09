import astropy.units as u
import astropy.coordinates as coord
from astropy.io import fits
import numpy as np
import pkg_resources
import os

class baumgardtMWGCs:
    """
    A class to handle data from the Baumgardt Milky Way Globular Clusters (MWGCs).

    This class provides methods to extract data from a FITS file, select a globular cluster of interest, 
    and build a covariance matrix of the observed kinematics of a given globular cluster.
    
    This path will store all the data, and return an astropy object of all the data in ICRS coordinates.

    Attributes
    ----------
    _pathtoclusterdata : str
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
            "rhopmrade":u.dimensionless_unscaled}
        ):
        """ extracts the data from the fits file and stores it in a dictionary

        """        
        
        # CHECK THAT THE FILE EXISTS
        self._pathtoclusterdata=pkg_resources.resource_filename('tstrippy', 'data/2023-03-28-merged.fits')
        self._unitkeys = ["RA","DEC","Rsun","RV","mualpha","mu_delta","Mass","rh_m","DM","rhopmrade"]
        self._unittypes= ["angle","angle","length","velocity","velocity","velocity","mass","length","length",""]
        
        assert self._pathtoclusterdata is not None, "The path to the cluster data is not valid"
        assert isinstance(self._pathtoclusterdata,str), "The path to the cluster data is not a string"
        if not os.path.exists(self._pathtoclusterdata):
            raise ValueError("The path to the cluster data is not valid")
        
        # CHECK THAT THE UNITS ARE VALID
        for key in units.keys():
            assert key in self._unitkeys, f"{key} is not a valid for the units dictionary"
        for key in self._unitkeys:
            assert key in units.keys(), f"{key} is not present in the units dictionary"

        # check that there are units present
        assert isinstance(units['RA'],u.UnitBase), f"{units[key]} is not a valid astropy unit"
        assert isinstance(units['DEC'],u.UnitBase), f"{units[key]} is not a valid astropy unit"
        assert isinstance(units['Rsun'],u.UnitBase), f"{units[key]} is not a valid astropy unit"
        assert isinstance(units['RV'],u.UnitBase), f"{units[key]} is not a valid astropy unit"
        assert isinstance(units['mualpha'],u.UnitBase), f"{units[key]} is not a valid astropy unit"
        assert isinstance(units['mu_delta'],u.UnitBase), f"{units[key]} is not a valid astropy unit"
        assert isinstance(units['Mass'],u.UnitBase), f"{units[key]} is not a valid astropy unit"
        assert isinstance(units['rh_m'],u.UnitBase), f"{units[key]} is not a valid astropy unit"
        assert isinstance(units['DM'],u.UnitBase), f"{units[key]} is not a valid astropy unit"
        
        # check that the are the correct units
        assert units['RA'].is_equivalent(u.degree), f"The units of RA is {units['RA']}, which is not equvalent to degree"
        assert units['DEC'].is_equivalent(u.degree), f"The units of DEC is {units['DEC']}, which is not equvalent to degree"
        assert units['Rsun'].is_equivalent(u.kpc), f"The units of Rsun is {units['Rsun']}, which is not equvalent to kpc"
        assert units['RV'].is_equivalent(u.km/u.s), f"The units of RV is {units['RV']}, which is not equvalent to km/s"
        assert units['mualpha'].is_equivalent(u.mas/u.yr), f"The units of mualpha is {units['mualpha']}, which is not equvalent to mas/yr"
        assert units['mu_delta'].is_equivalent(u.mas/u.yr), f"The units of mu_delta is {units['mu_delta']}, which is not equvalent to mas/yr"
        assert units['Mass'].is_equivalent(u.Msun), f"The units of Mass is {units['Mass']}, which is not equvalent to Msun"
        assert units['rh_m'].is_equivalent(u.kpc), f"The units of rh_m is {units['rh_m']}, which is not equvalent to kpc"
        assert units['DM'].is_equivalent(u.Msun), f"The units of DM is {units['DM']}, which is not equvalent to Msun"
        assert units['rhopmrade'].is_equivalent(u.dimensionless_unscaled), f"The units of rhopmrade is {units['rhopmrade']}, which is not equvalent to dimensionless_unscaled"

        self.units=units
        self._extractdatafromfits()


    def _extractdatafromfits(self):
        """ internal function to extract the data from the fits file and store it in a dictionary
        """        
        GCfits=fits.open(self._pathtoclusterdata)
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


