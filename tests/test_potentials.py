# tstrippy/tests/test_potentials.py
import tstrippy
import tstrippy.lib 
import tstrippy.dataParsers


def test_instantiation():
    print(tstrippy.lib.potentials)
    return None


def test_obtaining_parameters():
    MWparams = tstrippy.dataParsers.potential_parameters.pouliasis2017pii()
    assert len(MWparams) == 11
    
    return None

def test_pouliasis2017pii():
    MWparams = tstrippy.dataParsers.potential_parameters.pouliasis2017pii()
    x,y,z = 1,2,3
    out=tstrippy.lib.potentials.pouliasis2017pii(MWparams,x,y,z)
    assert(len(out)==4)
    return None