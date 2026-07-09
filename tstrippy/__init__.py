"""
TIDAL-STRIPPING-PYTHON
./tstrippy/__init__.py
"""
# starting from scratch 
from . import code
__all__ = ['code']

# from importlib import import_module
# import warnings

# # # Try to import Fortran modules (they're compiled into lib/)
# # # If they don't exist, provide helpful error messages
# # def _load_fortran_entry(module_name, attr_name):
# #     try:
# #         mod = import_module(f"{__name__}.lib.{module_name}")
# #         return getattr(mod, attr_name)
# #     except (ModuleNotFoundError, AttributeError):
# #         warnings.warn(
# #             f"Fortran module 'tstrippy.lib.{module_name}' not found. "
# #             "Have you built the package? Run: conda run -n tstrippy ./build.sh"
# #         )
# #         return None


# # simulator = _load_fortran_entry("simulator", "simulator")
# # gravity = _load_fortran_entry("gravity", "gravity")
# # mathutils = _load_fortran_entry("mathutils", "mathutils")

# # # try:
# # #     _gravity_ext = import_module(f"{__name__}.lib.gravity")
# # #     sphericalharmonicsbfe = _gravity_ext.sphericalharmonicsbfe
# # #     besselbfe = _gravity_ext.besselbfe
# # # except (ModuleNotFoundError, AttributeError):
# # #     sphericalharmonicsbfe = None
# # #     besselbfe = None
# # #     warnings.warn(
# # #         "Fortran backend modules are not available from 'tstrippy.lib.gravity'. "
# # #         "Have you built the package? Run: conda run -n tstrippy ./build.sh"
# # #     )

# # Import pure Python modules
# from . import io
# from . import code
# from . import simulator 
# from . import mathutils 
# from . import gravity

# # Define what's available at the top level
# __all__ = [
#     'simulator',
#     'gravity',
#     'mathutils',
#     'io',
#     'code',
# ]

# # Check for Fortran compiler
# import subprocess
# def _check_fortran_compiler():
#     try:
#         subprocess.run(['gfortran', '--version'], capture_output=True, check=False)  # noqa: F841
#     except FileNotFoundError:
#         warnings.warn(
#             "No Fortran compiler found. Some features of tstrippy may not work. "
#             "Please install gfortran 11+ for full functionality."
#         )

# _check_fortran_compiler()