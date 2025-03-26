project('tstrippy', 'fortran',
    version: '0.0.1',
    license: 'MIT'
)

# Define package name variable to avoid repetition
pkg_name = 'tstrippy'


py = import('python').find_installation()
fc_args = ['-fPIC']

# Get the Python site-packages directory
python_code = '''
import sys
from distutils.sysconfig import get_python_lib
print(get_python_lib())
'''

# Keep your existing install_dir setup for site-packages directory
install_dir = run_command('python', '-c', python_code, check: true).stdout().strip()
message('Python install dir: ' + install_dir)

# Define package name variable to avoid repetition
pkg_name = 'tstrippy'

# For extension modules, install them to the 'lib' subdirectory of your package
constants_ext = py.extension_module('constants',
    sources: 'tstrippy/src/constants.f90',
    install: true,
    install_dir: install_dir / pkg_name / 'lib',
    fortran_args: fc_args
)

potentials_ext = py.extension_module('potentials',
    sources: [
        'tstrippy/src/constants.f90', 
        'tstrippy/src/potentials.f90'
    ],
    install: true,
    install_dir: install_dir / pkg_name / 'lib',
    fortran_args: fc_args
)

integrator_ext = py.extension_module('integrator',
    sources: [
        'tstrippy/src/constants.f90',
        'tstrippy/src/potentials.f90',
        'tstrippy/src/perturbers.f90',
        'tstrippy/src/hostperturber.f90',
        'tstrippy/src/galacticbar.f90',
        'tstrippy/src/integrator.f90'
    ],
    install: true,
    install_dir: install_dir / pkg_name / 'lib',
    fortran_args: fc_args
)

# Then install Python files to their proper directories
py.install_sources(
    files(
        'tstrippy/__init__.py',
    ),
    subdir: install_dir / pkg_name
)

py.install_sources(
    files(
        'tstrippy/lib/__init__.py',
    ),
    subdir: install_dir / pkg_name / 'lib'
)

py.install_sources(
    files(
        'tstrippy/Parsers/__init__.py',
        'tstrippy/Parsers/baumgardtMWGCs.py',    
        'tstrippy/Parsers/potential_parameters.py',
    ),
    subdir: install_dir / pkg_name / 'Parsers'
)

py.install_sources(
    files(
        'tstrippy/code/__init__.py',
        'tstrippy/code/ergodic.py',
    ),
    subdir: install_dir/ pkg_name / 'code'
)

# And data files
install_data(
    files('tstrippy/data/2023-03-28-merged.fits',
          'tstrippy/data/MWrefframe001.yaml',
          'tstrippy/data/pouliasis2017pii.yaml',
          'tstrippy/data/unit_basis.yaml'),
    install_dir: install_dir / pkg_name / 'data'
)