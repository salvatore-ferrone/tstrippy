import tstrippy
import tstrippy.dataParsers

def test_pouliasis2017pii():
    myparams = tstrippy.dataParsers.potential_parameters.pouliasis2017pii()
    assert(len(myparams) == 11)