"""
tstrippy/lib/__init__.py

This directory contains compiled Fortran extension modules (simulator, gravity, mathutils).
They are imported in tstrippy/__init__.py and made available at the package top level.
"""
import warnings
from pkgutil import extend_path
from pathlib import Path
import sys


# Allow `tstrippy.lib` to resolve compiled extension modules that were installed
# into site-packages, even when importing from a source checkout.
__path__ = extend_path(__path__, __name__)

# When importing from a source checkout, compiled extensions may live in an
# installed site-packages copy of `tstrippy/lib`. Add those paths so
# `tstrippy.lib.simulator` etc. can still be resolved.
for _base in sys.path:
    try:
        _candidate = (Path(_base) / "tstrippy" / "lib").resolve()
    except Exception:
        continue
    if not _candidate.exists():
        continue
    _candidate_s = str(_candidate)
    if _candidate_s not in __path__:
        __path__.append(_candidate_s)

# Lazy loading: compiled modules are available as submodules when needed
# tstrippy/__init__.py handles the top-level imports with proper error handling

__all__ = []

def __getattr__(name):
    """Provide helpful error message if someone tries to access these directly."""
    if name in ['simulator', 'gravity', 'mathutils', "besselbfe", "sphericalharmonicsbfe"]:
        raise ModuleNotFoundError(
            f"Compiled Fortran module 'tstrippy.lib.{name}' not found. "
            f"Access it via 'tstrippy.{name}' instead (use the top-level import). "
            f"If not yet built, run: python -m pip install -e . --no-build-isolation"
        )
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")