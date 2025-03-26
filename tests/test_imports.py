import sys
import types

def test_top_level_modules():
    """Test that only the desired modules are visible at the top level."""
    import tstrippy
    
    # Check for expected modules
    assert hasattr(tstrippy, 'integrator'), "integrator module missing"
    assert hasattr(tstrippy, 'potentials'), "potentials module missing"
    assert hasattr(tstrippy, 'Parsers'), "Parsers module missing"
    assert hasattr(tstrippy, 'ergodic'), "ergodic module missing"
    
    
    # Check only expected modules are in top level (plus dunder methods)
    visible_modules = [x for x in dir(tstrippy) if not x.startswith('_')]
    expected_modules = ['integrator', 'potentials', 'Parsers', 'ergodic', 'lib', 'code']
    for module in expected_modules:
        assert module in visible_modules, f"{module} should be visible"
    
    # No unexpected modules
    unexpected = [x for x in visible_modules if x not in expected_modules]
    assert len(unexpected) == 0, f"Unexpected modules found: {unexpected}"

def test_module_types():
    """Test that the modules are of the correct types."""
    import tstrippy
    
    # Test Fortran modules
    assert 'leapfrogtofinalpositions' in dir(tstrippy.integrator), "integrator should expose integrator functions"
    assert 'plummer' in dir(tstrippy.potentials), "potentials should expose potentials functions"
    
    # Test ergodic is a proper module
    assert isinstance(tstrippy.ergodic, types.ModuleType), "ergodic should be a module"
    
    # Test ergodic contains expected functions
    ergodic_functions = [
        'DFPlummer', 'EscapeVelocity', 'PlummerMassProfile',
        'PlummerPotentialEnergy', 'PlummerRadius', 'UniformSphere',
        'convertHalfMassRadiusToPlummerRadius', 'isotropicplummer',
        'velocityCDF', 'velocitySampling'
    ]
    for func in ergodic_functions:
        assert hasattr(tstrippy.ergodic, func), f"ergodic.{func} function missing"

def test_direct_imports():
    """Test that modules can be directly imported."""
    # These should work without errors
    from tstrippy import integrator
    from tstrippy import potentials
    from tstrippy import ergodic
    from tstrippy import Parsers
    
    # Check that they're the same objects as accessed through tstrippy
    import tstrippy
    assert integrator is tstrippy.integrator
    assert potentials is tstrippy.potentials
    assert ergodic is tstrippy.ergodic
    assert Parsers is tstrippy.Parsers

def test_module_paths():
    """Test that modules are registered in sys.modules with the correct paths."""
    import tstrippy
    
    # Check expected modules are in sys.modules
    assert 'tstrippy' in sys.modules
    assert 'tstrippy.Parsers' in sys.modules
    assert 'tstrippy.code.ergodic' in sys.modules
    assert 'tstrippy.lib' in sys.modules
    

    