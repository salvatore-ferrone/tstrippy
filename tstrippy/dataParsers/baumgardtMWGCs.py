import astropy.units as u
from astropy.io import fits
import numpy as np
import pkg_resources


class baumgardtMWGCs:
    """Class for parsing the Baumgardt et al. 20XX-204ever catalog of Milky Way globular clusters
    
    NOTE. the catalog is available at https://people.smp.uq.edu.au/HolgerBaumgardt/globular/
    it is difficult to parse, I did this by hand and stored it in a fits file.
    
    
    TO DO LIST:
        Generate a code that queries the catalog, downloads the fits file and parses it.
        - sample from skewed distribution of the distances
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

        :param units: an optional dictionary with the requested units for each field. The keys must be correct and the values must be astropy units that are compatible with those in the fits file, defaults to { "RA":u.degree, "DEC":u.degree, "Rsun":u.kpc, "RV":u.km/u.s, "mualpha":u.mas/u.yr, "mu_delta":u.mas/u.yr, "Mass":u.Msun, "rh_m":u.kpc, "DM":u.Msun, "rhopmrade":""}
        :type units: dict, optional
        """        
        
        self.pathtoclusterdata=pkg_resources.resource_filename('tstrippy', 'data/2023-03-28-merged.fits')
        
        self._extractdatafromfits()
        self.units=units


    def _extractdatafromfits(self):
        """
        parses the fits file and stores the data in a dictionary
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
        """
        Make the properties of the GC of interest available for the class
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
        """
        Purpose:
            Build the co-variance matrix of the observed kinematics of a given GC
        Input:
            dataTable: the astropy table of the GC data
            GCidx: the index of the GC in the table
        Output:
            means: the mean values of the kinematics [RA,DEC,Rsun,RV,mualpha,mu_delta]
            covMat: the co-variance matrix

        NOTE: the order of the kinematics is:
            [RA,DEC,Rsun,RV,mualpha,mu_delta]
        NOTE:
            the errors on RA and DEC are set to zero
            only the covariance between the proper motions is given
            intendned to be sampled with 
        EXAMPLE:
        self.BuildCovarianceMatrix()
        RA,DEC,D,RV,mualpha,mu_delta,Mass,halfmassradius = np.random.multivariate_normal(self.meanParameters,self.covarianceParameters)
        """
        # 
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


