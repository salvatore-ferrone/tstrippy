import tstrippy
import numpy as np
from astropy import units as u
# test initialization
def test_instantiation():
    myobj = tstrippy.Parsers.baumgardtMWGCs()
    assert isinstance(myobj, tstrippy.Parsers.baumgardtMWGCs)
    return None

def test_initialization():
    myobj = tstrippy.Parsers.baumgardtMWGCs()
    assert isinstance(myobj.units, dict)
    assert isinstance(myobj._unitkeys, list)
    for key in myobj._unitkeys:
        assert key in myobj.units.keys(), f"{key} is in myobj._unitkeys but not in myobj.units"
    for key in myobj.units.keys():
        assert key in myobj._unitkeys, f"{key} is in myobj.units but not in myobj._unitkeys"
    assert isinstance(myobj.units['RA'], u.UnitBase)
    assert isinstance(myobj.units['DEC'], u.UnitBase)
    assert isinstance(myobj.units['Rsun'], u.UnitBase)
    assert isinstance(myobj.units['RV'], u.UnitBase)
    assert isinstance(myobj.units['mualpha'], u.UnitBase)
    assert isinstance(myobj.units['mu_delta'], u.UnitBase)
    assert isinstance(myobj.units['Mass'], u.UnitBase)
    assert isinstance(myobj.units['rh_m'], u.UnitBase)
    assert isinstance(myobj.units['rhopmrade'], u.UnitBase)
    assert isinstance(myobj._pathtoclusterdata, str)
    assert isinstance(myobj._fitskeys, list)
    assert isinstance(myobj.data, dict)
    for key in myobj._fitskeys:
        assert key in myobj.data.keys(), f"{key} is in myobj._fitskeys but not in myobj.data"
    for key in myobj.data.keys():
        assert key in myobj._fitskeys, f"{key} is in myobj.data but not in myobj._fitskeys"
    
    return None


def test_getGCCovarianceMatrix():
    ''' test that the covariance matrix is correct by comparing the numerical mean and standard deviation to the true mean and standard deviation
    the difference should be less than 5/sqrt(N) where N is the number of samples.
    i.e. I am comparing the difference of the numerical mean and the true mean to the standard error of the mean.
    if the standard deviation is zero, then the difference should be less than 1e-11, which should capture numerical error.'''
    numericalerror=1e-11
    GCnames=["Pal5","NGC3201", "Pal10", "NGC5466", "NGC4147", "NGC4590", "NGC5024", "NGC5053", "NGC5272"]
    nsamps=int(1e2)
    tolerance=5
    myobj = tstrippy.Parsers.baumgardtMWGCs()
    # pick a random GC name
    
    for GCname in GCnames:
        means,cov=myobj.getGCCovarianceMatrix(GCname)
        assert isinstance(means, np.ndarray)
        assert isinstance(cov, np.ndarray)
        assert cov.shape == (8,8), f"cov.shape={cov.shape} but should be (8,8)"

        
        for column in range(cov.shape[0]):
            assertstring="\n\tGCname={GCname}\n\tcolumn={column}\n\tnsamps={nsamps}".format(GCname=GCname,column=column,nsamps=nsamps)
            samples=np.random.multivariate_normal(means,cov,nsamps)
            trueMean=means[column]
            trueSTD=np.sqrt(cov[column,column])
            numericalSTD = np.std(samples[:,column])
            numericalMean=np.mean(samples[:,column])
            meandiff=np.abs(numericalMean-trueMean)
            additionalinfo="\n\ttrueMean={trueMean}\n\ttrueSTD={trueSTD}\n\tnumericalMean={numericalMean}\n\tnumericalSTD={numericalSTD}".format(trueMean=trueMean,numericalMean=numericalMean,trueSTD=trueSTD,numericalSTD=numericalSTD)
            if trueSTD==0:
                assert np.abs(numericalMean-trueMean) < numericalerror, "trueSTD==0 but np.abs(numericalMean-trueMean)={}".format(np.abs(numericalMean-trueMean)) + assertstring + additionalinfo
            else:
                assert meandiff/numericalSTD < tolerance/np.sqrt(nsamps), "meandiff/numericalSTD={} < {}".format(meandiff/numericalSTD,tolerance/np.sqrt(nsamps)) + assertstring+additionalinfo
                    


    return None