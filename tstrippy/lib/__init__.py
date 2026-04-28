"""
tstrippy/lib/__init__.py

This directory contains compiled Fortran extension modules (integrator, potentials, mathutils).
They are imported in tstrippy/__init__.py and made available at the package top level.
"""
import warnings

# Lazy loading: compiled modules are available as submodules when needed
# tstrippy/__init__.py handles the top-level imports with proper error handling

__all__ = []

def __getattr__(name):
    """Provide helpful error message if someone tries to access these directly."""
    if name in ['integrator', 'potentials', 'mathutils']:
        raise ModuleNotFoundError(
            f"Compiled Fortran module 'tstrippy.lib.{name}' not found. "
            f"Access it via 'tstrippy.{name}' instead (use the top-level import). "
            f"If not yet built, run: python -m pip install -e . --no-build-isolation"
        )
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")