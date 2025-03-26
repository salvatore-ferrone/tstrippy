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