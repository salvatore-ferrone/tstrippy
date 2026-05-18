"""tstrippy/io/__init__.py"""
import warnings

try:
    from . import mwgcs
except Exception as exc:
    mwgcs = None
    warnings.warn(f"Failed to import tstrippy.io.mwgcs: {exc}")
try:
    from . import milkyway_models
except Exception as exc:
    milkyway_models = None
    warnings.warn(f"Failed to import tstrippy.io.milkyway_models: {exc}")
try:
    from . import reference_frames
except Exception as exc:
    reference_frames = None
    warnings.warn(f"Failed to import tstrippy.io.reference_frames: {exc}")
from .write_simulation_hdf5 import write_simulation_hdf5, read_snapshot_binary, read_orbit_binary

__all__ = ["mwgcs", "milkyway_models", "reference_frames",
           "write_simulation_hdf5", "read_snapshot_binary", "read_orbit_binary"]