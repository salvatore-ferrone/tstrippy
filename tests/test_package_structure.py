import importlib
import sys
import types


def test_top_level_public_api():
    import tstrippy

    expected = {
        "simulator",
        "gravity",
        "mathutils",
        "besselbfe",
        "sphericalharmonicsbfe",
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

    assert hasattr(tstrippy.simulator, "setinitialconditions")
    assert hasattr(tstrippy.simulator, "setscheme")
    assert hasattr(tstrippy.simulator, "settimestamps")
    assert hasattr(tstrippy.simulator, "setbackwardorbit")
    assert hasattr(tstrippy.simulator, "trim_orbits")
    assert hasattr(tstrippy.simulator, "clear")
    assert hasattr(tstrippy.simulator, "cleargravitycomponents")
    assert hasattr(tstrippy.simulator, "set_gravitational_constant")
    assert hasattr(tstrippy.simulator, "add_component")
    assert hasattr(tstrippy.simulator, "finalize")
    assert hasattr(tstrippy.simulator, "run")
    assert hasattr(tstrippy.simulator, "initwritesnapshots")
    assert hasattr(tstrippy.simulator, "initwriteorbits")
    assert hasattr(tstrippy.simulator, "write_snapshot_file")
    assert hasattr(tstrippy.simulator, "open_orbit_files")
    assert hasattr(tstrippy.simulator, "close_orbit_files")
    assert hasattr(tstrippy.simulator, "write_orbit_records")
    assert hasattr(tstrippy.simulator, "build_fixed_timestamps")
    assert hasattr(tstrippy.simulator, "compute_memory_particle_limit")
    assert hasattr(tstrippy.simulator, "allocate_orbits")
    assert hasattr(tstrippy.simulator, "leapfrog")
    assert hasattr(tstrippy.simulator, "forest_ruth")



def test_code_namespace_structure():
    import tstrippy

    assert set(tstrippy.code.__all__) == {"bfe", "sampling", "orbits"}


def test_io_namespace_structure():
    import tstrippy

    assert isinstance(tstrippy.io, types.ModuleType)
    assert set(tstrippy.io.__all__) == {
        "mwgcs",
        "milkyway_models",
        "reference_frames",
        "write_simulation_hdf5",
    }

    assert isinstance(tstrippy.io.mwgcs, types.ModuleType)
    assert isinstance(tstrippy.io.milkyway_models, types.FunctionType)
    assert isinstance(tstrippy.io.reference_frames, types.FunctionType)
    assert isinstance(tstrippy.io.write_simulation_hdf5, types.FunctionType)


def test_direct_import_paths():
    modules = {
        "tstrippy",
        "tstrippy.io",
        "tstrippy.io.mwgcs",
        "tstrippy.io.milkyway_models",
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
    _ = tstrippy.io.mwgcs
    _ = tstrippy.io.milkyway_models
    _ = tstrippy.io.write_simulation_hdf5
    _ = tstrippy.io.reference_frames

    expected = {
        "tstrippy",
        "tstrippy.lib",
        "tstrippy.io",
        "tstrippy.io.mwgcs",
        "tstrippy.io.milkyway_models",
        "tstrippy.io.reference_frames",
        "tstrippy.io.write_simulation_hdf5",
        "tstrippy.code",
        "tstrippy.code.bfe",
        "tstrippy.code.sampling",
        "tstrippy.code.orbits",
    }

    assert expected.issubset(set(sys.modules))