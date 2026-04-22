import importlib
import sys
import types


def test_top_level_public_api():
    import tstrippy

    expected = {
        "integrator",
        "potentials",
        "mathutils",
        "io",
        "code",
    }

    for name in expected:
        assert hasattr(tstrippy, name), f"missing top-level attribute: {name}"

    visible = {name for name in dir(tstrippy) if not name.startswith("_")}
    assert expected.issubset(visible)


def test_top_level_module_types():
    import tstrippy

    assert isinstance(tstrippy.io, types.ModuleType)
    assert isinstance(tstrippy.code, types.ModuleType)
    assert isinstance(tstrippy.code.bfe, types.ModuleType)
    assert isinstance(tstrippy.code.sampling, types.ModuleType)
    assert isinstance(tstrippy.code.orbits, types.ModuleType)

    assert hasattr(tstrippy.integrator, "leapfrogtofinalpositions")
    assert hasattr(tstrippy.potentials, "plummer")
    assert hasattr(tstrippy.mathutils, "linear_interp_scalar")


def test_code_namespace_structure():
    import tstrippy

    assert set(tstrippy.code.__all__) == {"bfe", "sampling", "orbits"}


def test_io_namespace_structure():
    import tstrippy

    assert isinstance(tstrippy.io, types.ModuleType)
    assert set(tstrippy.io.__all__) == {
        "baumgardt_gcs",
        "potential_parameters",
        "reference_frames",
    }

    assert isinstance(tstrippy.io.baumgardt_gcs, types.ModuleType)
    assert isinstance(tstrippy.io.potential_parameters, types.ModuleType)
    assert isinstance(tstrippy.io.reference_frames, types.ModuleType)


def test_direct_import_paths():
    modules = {
        "tstrippy",
        "tstrippy.io",
        "tstrippy.io.baumgardt_gcs",
        "tstrippy.io.potential_parameters",
        "tstrippy.io.reference_frames",
        "tstrippy.code",
        "tstrippy.code.bfe",
        "tstrippy.code.sampling",
        "tstrippy.code.orbits",
    }

    for module_name in modules:
        module = importlib.import_module(module_name)
        assert isinstance(module, types.ModuleType), f"failed import: {module_name}"


def test_sys_modules_registration():
    import tstrippy

    _ = tstrippy.io
    _ = tstrippy.code
    _ = tstrippy.code.bfe
    _ = tstrippy.code.sampling
    _ = tstrippy.code.orbits
    _ = tstrippy.io.baumgardt_gcs
    _ = tstrippy.io.potential_parameters
    _ = tstrippy.io.reference_frames

    expected = {
        "tstrippy",
        "tstrippy.lib",
        "tstrippy.io",
        "tstrippy.io.baumgardt_gcs",
        "tstrippy.io.potential_parameters",
        "tstrippy.io.reference_frames",
        "tstrippy.code",
        "tstrippy.code.bfe",
        "tstrippy.code.sampling",
        "tstrippy.code.orbits",
    }

    assert expected.issubset(set(sys.modules))