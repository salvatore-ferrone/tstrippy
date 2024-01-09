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
    data : dict
        The data extracted from the FITS file.
    units : dict
        The units for the data fields.
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
            "rhopmrade":u.dimensionless_unscaled}
        ):
        """ extracts the data from the fits file and stores it in a dictionary

        """        
        
        # CHECK THAT THE FILE EXISTS
        self._pathtoclusterdata=pkg_resources.resource_filename('tstrippy', 'data/2023-03-28-merged.fits')
        self._unitkeys = ["RA","DEC","Rsun","RV","mualpha","mu_delta","Mass","rh_m","rhopmrade"]
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
        
        # check that the are the correct units
        assert units['RA'].is_equivalent(u.degree), f"The units of RA is {units['RA']}, which is not equvalent to degree"
        assert units['DEC'].is_equivalent(u.degree), f"The units of DEC is {units['DEC']}, which is not equvalent to degree"
        assert units['Rsun'].is_equivalent(u.kpc), f"The units of Rsun is {units['Rsun']}, which is not equvalent to kpc"
        assert units['RV'].is_equivalent(u.km/u.s), f"The units of RV is {units['RV']}, which is not equvalent to km/s"
        assert units['mualpha'].is_equivalent(u.mas/u.yr), f"The units of mualpha is {units['mualpha']}, which is not equvalent to mas/yr"
        assert units['mu_delta'].is_equivalent(u.mas/u.yr), f"The units of mu_delta is {units['mu_delta']}, which is not equvalent to mas/yr"
        assert units['Mass'].is_equivalent(u.Msun), f"The units of Mass is {units['Mass']}, which is not equvalent to Msun"
        assert units['rh_m'].is_equivalent(u.kpc), f"The units of rh_m is {units['rh_m']}, which is not equvalent to kpc"
        assert units['rhopmrade'].is_equivalent(u.dimensionless_unscaled), f"The units of rhopmrade is {units['rhopmrade']}, which is not equvalent to dimensionless_unscaled"
        
        self.units=units
        self._extractdatafromfits()


    def _extractdatafromfits(self):
        """ internal function to extract the data from the fits file and store it in a dictionary
        """  
        
        
        try:
            GCfits = fits.open(self._pathtoclusterdata,'readonly')
        except FileNotFoundError:
            raise FileNotFoundError(f"The file {self._pathtoclusterdata} does not exist.")
        except OSError:
            raise OSError(f"The file {self._pathtoclusterdata} is not a valid FITS file.")
        
        
        num_tunits = sum(1 for key in GCfits[1].header if key.startswith('TUNIT'))
        data = {}
        for n in range(1,num_tunits+1):
            myttype=GCfits[1].header['TTYPE'+str(n)]
            mytunit=GCfits[1].header['TUNIT'+str(n)]
            if mytunit=='Name':
                data[myttype]=GCfits[1].data[GCfits[1].header['TTYPE'+str(n)]]
            else:
                data[myttype]=GCfits[1].data[GCfits[1].header['TTYPE'+str(n)]]*u.Unit(GCfits[1].header['TUNIT'+str(n)])

        
        # convert to the desired units
        data['RA']=data['RA'].to(self.units['RA'])
        data['DEC']=data['DEC'].to(self.units['DEC'])
        data['Rsun']=data['Rsun'].to(self.units['Rsun'])
        data['ERsun']=data['ERsun'].to(self.units['Rsun'])
        data['RV']=data['RV'].to(self.units['RV'])
        data['ERV']=data['ERV'].to(self.units['RV'])
        data['mualpha']=data['mualpha'].to(self.units['mualpha'])
        data['ERmualpha']=data['ERmualpha'].to(self.units['mualpha'])
        data['mu_delta']=data['mu_delta'].to(self.units['mu_delta'])
        data['ERmu_delta']=data['ERmu_delta'].to(self.units['mu_delta'])
        data['Mass']=data['Mass'].to(self.units['Mass'])
        data['DM']=data['DM'].to(self.units['Mass'])
        data['rh_m']=data['rh_m'].to(self.units['rh_m'])
        data['rhopmrade']=data['rhopmrade'].to(self.units['rhopmrade'])
        
        self.data=data
        
        GCfits.close()
    
    
    def getGCCovarianceMatrix(self,GCname):
        """Pick a globular cluster of interest and extract its data to tbe class directly

        Parameters
        ----------
        GCname : str
            name of the globular cluter of interest. Must be one of the names in the fits file
        
        Returns
        -------
        means : list
            list of the means of the parameters
        covarianceMatrix : np.array
            covariance matrix of the parameters
        
        Note
        ----
        the order of the data are:
            [RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m]
        Example
        --------
        means,covarianceMatrix = GC.getGCCovarianceMatrix('NGC104')
        
        """        
        assert isinstance(GCname,str), f"{GCname} is not a string"
        assert GCname in self.data['Cluster'], f"{GCname} is not a valid globular cluster name"
        
        GCdex=np.where(self.data['Cluster']==GCname)[0][0]
        
        # drop the units for numpy
        RA=self.data['RA'][GCdex].value
        ERRA=0 # no error
        DEC=self.data['DEC'][GCdex].value
        ERDEC=0 # no error
        Rsun=self.data['Rsun'][GCdex].value
        ERsun=self.data['ERsun'][GCdex].value
        RV=self.data['RV'][GCdex].value
        ERV=self.data['ERV'][GCdex].value
        mualpha=self.data['mualpha'][GCdex].value
        ERmualpha=self.data['ERmualpha'][GCdex].value
        mu_delta=self.data['mu_delta'][GCdex].value
        ERmu_delta=self.data['ERmu_delta'][GCdex].value
        rhopmrade=self.data['rhopmrade'][GCdex].value
        Mass=self.data['Mass'][GCdex].value
        DM=self.data['DM'][GCdex].value
        rh_m=self.data['rh_m'][GCdex].value
        ERrh_m=0 # no error
        
        means=[RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m]
        sigmas=[ERRA,ERDEC,ERsun,ERV,ERmualpha,ERmu_delta,DM,ERrh_m]
        Nparams = len(means)
        # put in the variances
        covarianceMatrix=np.zeros((Nparams,Nparams))
        for i in range(Nparams):
            covarianceMatrix[i,i]=sigmas[i]**2
        # convert correlation to covariance and put in the covariance
        covarianceMatrix[4,5]=covarianceMatrix[5,4]=rhopmrade*ERmualpha*ERmu_delta # only the covariance between the proper motions
        return means,covarianceMatrix
        

