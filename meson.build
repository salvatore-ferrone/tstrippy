project('tstrippy', 'fortran',
    version: '0.0.1',
    license: 'MIT'
)

# Define package name variable to avoid repetition
pkg_name = 'tstrippy'

# Find Python and f2py
py = import('python').find_installation()
f2py = find_program('f2py3', 'f2py')
fc_args = ['-fPIC']

# Get the Python site-packages directory
python_code = '''
import sysconfig
print(sysconfig.get_path('purelib'))
'''

# Get Python extension suffix for module naming
ext_suffix_code = '''
import sysconfig
print(sysconfig.get_config_var('EXT_SUFFIX'))
'''

install_dir = run_command(py, '-c', python_code, check: true).stdout().strip()
ext_suffix = run_command(py, '-c', ext_suffix_code, check: true).stdout().strip()
message('Python install dir: ' + install_dir)
message('Python extension suffix: ' + ext_suffix)




message('Running direct meson build')
# Paths for direct meson install (absolute paths)
module_install_dir = install_dir / pkg_name / 'lib'
python_install_subdir = install_dir / pkg_name


# Setup environment for f2py
fc_env = environment()
fc_env.set('FC', 'gfortran')
fc_env.set('FFLAGS', '-Wno-unused-variable -fPIC')


constants_sources = ['tstrippy/src/constants.f90']
constants_out = 'constants' + ext_suffix
constants_build = custom_target('constants',
    output: constants_out,
    input: constants_sources,
    command: [f2py, '-c', '@INPUT@', '--backend=meson', '--quiet', '-m', 'constants'],
    env: fc_env,
    install: true,
    install_dir: module_install_dir
)

potentials_sources = [
    'tstrippy/src/constants.f90',
    'tstrippy/src/potentials.f90'
]
potentials_out = 'potentials' + ext_suffix

potentials_build = custom_target('potentials',
    output: potentials_out,
    input: potentials_sources,
    command: [f2py, '-c', '@INPUT@', '--backend=meson', '--quiet', '-m', 'potentials'],
    env: fc_env,
    install: true,
    install_dir: module_install_dir
)

integrator_sources = [
    'tstrippy/src/constants.f90',
    'tstrippy/src/potentials.f90',
    'tstrippy/src/perturbers.f90',
    'tstrippy/src/hostperturber.f90',
    'tstrippy/src/galacticbar.f90',
    'tstrippy/src/integrator.f90'
]
integrator_out = 'integrator' + ext_suffix
integrator_build = custom_target('integrator',
    output: integrator_out,
    input: integrator_sources,
    command: [f2py, '-c', '@INPUT@', '--backend=meson', '--quiet', '-m', 'integrator'],
    env: fc_env,
    install: true,
    install_dir: module_install_dir
)

# Then install Python files to their proper directories
py.install_sources(
    files(
        'tstrippy/__init__.py',
    ),
    subdir: python_install_subdir
)

py.install_sources(
    files(
        'tstrippy/lib/__init__.py',
    ),
    subdir: python_install_subdir / 'lib'
)

py.install_sources(
    files(
        'tstrippy/Parsers/__init__.py',
        'tstrippy/Parsers/baumgardtMWGCs.py',    
        'tstrippy/Parsers/potential_parameters.py',
    ),
    subdir: python_install_subdir / 'Parsers'
)

py.install_sources(
    files(
        'tstrippy/code/__init__.py',
        'tstrippy/code/ergodic.py',
    ),
    subdir: python_install_subdir / 'code'
)

# And data files
install_data(
    files('tstrippy/data/2023-03-28-merged.fits',
          'tstrippy/data/MWrefframe001.yaml',
          'tstrippy/data/pouliasis2017pii.yaml',
          'tstrippy/data/unit_basis.yaml'),
    install_dir: python_install_subdir / 'data'
)