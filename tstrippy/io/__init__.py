"""tstrippy/io/__init__.py"""
import warnings

try:
    from . import mwgcs
except Exception as exc:
    mwgcs = None
    warnings.warn(f"Failed to import tstrippy.io.mwgcs: {exc}")
try:
    from . import milkyway_models as _milkyway_models_module
    milkyway_models = _milkyway_models_module.milkyway_models
    # Attach helper methods to keep a compact callable UX.
    milkyway_models.get_model = _milkyway_models_module.get_model
    milkyway_models.save_model_yaml = _milkyway_models_module.save_model_yaml
    milkyway_models._to_yaml_compatible = _milkyway_models_module._to_yaml_compatible
except Exception as exc:
    milkyway_models = None
    warnings.warn(f"Failed to import tstrippy.io.milkyway_models: {exc}")
try:
    from . import reference_frames as _reference_frames_module
    reference_frames = _reference_frames_module.reference_frames
    # Attach helper methods to keep a compact callable UX.
    reference_frames.get_reference_frame = _reference_frames_module.get_reference_frame
    reference_frames.available_reference_frames = _reference_frames_module.available_reference_frames
    reference_frames.load_yaml = _reference_frames_module.load_yaml
    reference_frames.MWrefframeFerrone2023 = _reference_frames_module.MWrefframeFerrone2023
except Exception as exc:
    reference_frames = None
    warnings.warn(f"Failed to import tstrippy.io.reference_frames: {exc}")
try:
    from . import write_simulation_hdf5 as _write_simulation_hdf5
    write_simulation_hdf5 = _write_simulation_hdf5.write_simulation_hdf5
    write_simulation_hdf5.read_snapshot_binary = _write_simulation_hdf5.read_snapshot_binary
    write_simulation_hdf5.read_orbit_binary = _write_simulation_hdf5.read_orbit_binary
except Exception as exc:
    write_simulation_hdf5 = None
    warnings.warn(f"Failed to import tstrippy.io.write_simulation_hdf5: {exc}")

from .write_simulation_hdf5 import write_simulation_hdf5, read_snapshot_binary, read_orbit_binary

__all__ = ["mwgcs", "milkyway_models", "reference_frames", "write_simulation_hdf5",]