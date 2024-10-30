import tstrippy
from astropy import units as u
from astropy import constants as const


def test_loadunits():
    # Load the units
    unitbasis = tstrippy.Parsers.potential_parameters.unitbasis
    unitT=u.Unit(unitbasis['time'])
    unitV=u.Unit(unitbasis['velocity'])
    unitD=u.Unit(unitbasis['distance'])
    unitM=u.Unit(unitbasis['mass'])
    unitG=u.Unit(unitbasis['G'])
    G = const.G.to(unitG).value

    assert isinstance(unitT, u.UnitBase)
    assert isinstance(unitV, u.UnitBase)
    assert isinstance(unitD, u.UnitBase)
    assert isinstance(unitM, u.UnitBase)
    assert isinstance(unitG, u.UnitBase)
    assert isinstance(G, float)

    return None
