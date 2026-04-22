"""
TIDAL-STRIPPING-PYTHON
./tstrippy/__init__.py
"""
from importlib import import_module

# Import the Fortran modules and move them to the top level
from .lib.integrator import integrator
from .lib.potentials import potentials
from .lib.mathutils import mathutils

# Import Parsers module
from . import io

# Define what's available at the top level
__all__ = [
    'integrator',
    'potentials',
    'mathutils',
    'io',
    'sampling',
    'bfe',
    "orbits"
]

# Check for Fortran compiler
import subprocess
import warnings
def _check_fortran_compiler():
    try:
        subprocess.run(['gfortran', '--version'], capture_output=True)
    except FileNotFoundError:
        warnings.warn(
            "No Fortran compiler found. Some features of tstrippy may not work. "
            "Please install gfortran 11+ for full functionality."
        )

_check_fortran_compiler()


def __getattr__(name):
    if name == 'bfe':
        return import_module('.code.bfe', __name__)
    if name == 'sampling':
        return import_module('.code.sampling', __name__)
    if name == 'orbits':
        return import_module('.code.orbits', __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

# delete subprocess and warnings
del subprocess, warnings 