project('tstrippy', 'fortran',
    version: '0.0.1',
    license: 'MIT'
)

py = import('python').find_installation()
fc_args = ['-fPIC']


#### NOTE 
# To ensure that the package is installed to the conda environment, I need to run this python code to get the install directory
# the default the comes with meson, `install_dir = py.get_install_dir()` will install the package to the system python environment
# which is not what I want and will not work because I don't have sudo access to the system python environment
python_code = '''
import sys
from distutils.sysconfig import get_python_lib
print(get_python_lib())
'''
install_dir = run_command('python', '-c', python_code, check: true).stdout().strip()

message('Python install dir: ' + install_dir)


# ---- Build Fortran Shared Libraries ----
constants_mod = shared_library('constants',
    sources: 'tstrippy/src/constants.f90',
    install: true,
    install_dir: install_dir / 'tstrippy/lib',
    fortran_args: fc_args
)

potentials_mod = shared_library('potentials',
    sources: 'tstrippy/src/potentials.f90',
    install: true,
    install_dir: install_dir / 'tstrippy/lib',
    fortran_args: fc_args
)

integrator_mod = shared_library('integrator',
    sources: [
        'tstrippy/src/constants.f90',
        'tstrippy/src/potentials.f90',
        'tstrippy/src/perturbers.f90',
        'tstrippy/src/hostperturber.f90',
        'tstrippy/src/galacticbar.f90',
        'tstrippy/src/integrator.f90'
    ],
    install: true,
    install_dir: install_dir / 'tstrippy/lib',
    fortran_args: fc_args
)

# ---- Install Python Package ----
# py.install_sources(meson.source_root() / 'tstrippy', subdir: 'tstrippy')
py.install_sources(
    files(
        'tstrippy/Parsers/baumgardtMWGCs.py',    
        'tstrippy/Parsers/potential_parameters.py',
        'tstrippy/Parsers/__init__.py',
        'tstrippy/lib/__init__.py',
        'tstrippy/__init__.py',
        'tstrippy/code/ergodic.py',
        'tstrippy/code/__init__.py'), 
    subdir: install_dir / 'tstrippy',
)

# ---- Install Data Files ----
install_data(
  files('tstrippy/data/2023-03-28-merged.fits',
        'tstrippy/data/MWrefframe001.yaml',
        'tstrippy/data/pouliasis2017pii.yaml',
        'tstrippy/data/unit_basis.yaml'),
    install_dir:install_dir / 'tstrippy/data'
)
# install_data(meson.source_root() / 'tstrippy/data/', install_dir: py.get_install_dir() / 'tstrippy/data')

# ---- Install Using pip (Optional) ----
# meson.add_install_script(py, '-m', 'pip', 'install', '--prefix', get_option('prefix'), '.')
