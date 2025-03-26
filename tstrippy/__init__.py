"""
TIDAL-STRIPPING-PYTHON
./tstrippy/__init__.py
"""
# Import the Fortran modules and move them to the top level
from .lib.integrator import integrator
from .lib.potentials import potentials

# Import the ergodic module and make it available at the top level
from .code import ergodic

# Import Parsers module
from . import Parsers

# Define what's available at the top level
__all__ = [
    'integrator',
    'potentials',
    'Parsers',
    'ergodic',
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