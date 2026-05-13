"""
TIDAL-STRIPPING-PYTHON
./tstrippy/__init__.py
"""
from importlib import import_module
import warnings

# Try to import Fortran modules (they're compiled into lib/)
# If they don't exist, provide helpful error messages
try:
    from .lib.simulator import simulator
except ModuleNotFoundError:
    simulator = None
    warnings.warn(
        "Fortran module 'tstrippy.lib.simulator' not found. "
        "Have you built the package? Run: python -m pip install -e . --no-build-isolation"
    )

try:
    from .lib.gravity import gravity
except ModuleNotFoundError:
    gravity = None
    warnings.warn(
        "Fortran module 'tstrippy.lib.gravity' not found. "
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

try:
    _gravity_ext = import_module(f"{__name__}.lib.gravity")
    sphericalharmonicsbfe = _gravity_ext.sphericalharmonicsbfe
    besselbfe = _gravity_ext.besselbfe
except (ModuleNotFoundError, AttributeError):
    sphericalharmonicsbfe = None
    besselbfe = None
    warnings.warn(
        "Fortran backend modules are not available from 'tstrippy.lib.gravity'. "
        "Have you built the package? Run: python -m pip install -e . --no-build-isolation"
    )

# Import pure Python modules
from . import io
from . import code

# Define what's available at the top level
__all__ = [
    'simulator',
    'gravity',
    'mathutils',
    "sphericalharmonicsbfe",
    "besselbfe",
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