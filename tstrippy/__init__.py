"""
TIDAL-STRIPPING-PYTHON
./tstrippy/__init__.py
"""
from importlib import import_module
import warnings

# Try to import Fortran modules (they're compiled into lib/)
# If they don't exist, provide helpful error messages
try:
    from .lib.integrator import integrator
except ModuleNotFoundError:
    integrator = None
    warnings.warn(
        "Fortran module 'tstrippy.lib.integrator' not found. "
        "Have you built the package? Run: python -m pip install -e . --no-build-isolation"
    )

try:
    from .lib.potentials import potentials
except ModuleNotFoundError:
    potentials = None
    warnings.warn(
        "Fortran module 'tstrippy.lib.potentials' not found. "
        "Have you built the package? Run: python -m pip install -e . --no-build-isolation"
    )

try:
    from .lib.mathutils import mathutils
except ModuleNotFoundError:
    mathutils = None
    warnings.warn(
        "Fortran module 'tstrippy.lib.mathutils' not found. "
        "Have you built the package? Run: python -m pip install -e . --no-build-isolation"
    )

# Import pure Python modules
from . import io
from . import code

# Define what's available at the top level
__all__ = [
    'integrator',
    'potentials',
    'mathutils',
    'io',
    'code',
]

# Check for Fortran compiler
import subprocess
def _check_fortran_compiler():
    try:
        subprocess.run(['gfortran', '--version'], capture_output=True, check=False)  # noqa: F841
    except FileNotFoundError:
        warnings.warn(
            "No Fortran compiler found. Some features of tstrippy may not work. "
            "Please install gfortran 11+ for full functionality."
        )

_check_fortran_compiler()