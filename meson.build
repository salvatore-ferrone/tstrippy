project('tstrippy', ['c'],
    version: '0.0.1',
    license: 'MIT'
)

# Define package name variable to avoid repetition
pkg_name = 'tstrippy'

# Find Python and f2py
py = import('python').find_installation()
# f2py = find_program('f2py') # could be system wide and incorrect
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

# Use Python to locate the correct f2py from NumPy
f2py_code = '''
import sys, os, numpy
f2py_path = os.path.join(os.path.dirname(numpy.__file__), 'f2py', 'f2py.py')
if os.path.exists(f2py_path):
    print('python:' + f2py_path)
else:
    import shutil
    f2py_bin = shutil.which('f2py')
    print(f2py_bin if f2py_bin else 'f2py')
'''

f2py_result = run_command(py, '-c', f2py_code, check: true).stdout().strip()
if f2py_result.startswith('python:')
    # Use Python to run f2py.py script
    f2py_script = f2py_result.split(':', 1)[1]
    f2py = [py, f2py_script]
else
    # Use binary f2py from PATH
    f2py = find_program(f2py_result)
endif

message('Running direct meson build')
# Paths for direct meson install (absolute paths)
module_install_dir = install_dir / pkg_name / 'lib'
python_install_subdir = install_dir / pkg_name


# Setup environment for f2py
fc_env = environment()
fc_env.set('FC', 'gfortran')
fc_env.set('FFLAGS', '-Wno-unused-variable -fPIC')
fc_env.set('CFLAGS', '-Wno-misleading-indentation')  # Suppress misleading indentation warning

# Suppress compact unwind warnings from macOS linker
if host_machine.system() == 'darwin'
    # For macOS, add linker flags to suppress unwind warnings
    fc_env.set('LDFLAGS', '-Wl,-no_compact_unwind')
endif

mathutils_sources = [
    'tstrippy/src/mathutils.f90',
]
mathutils_out = 'mathutils' + ext_suffix
mathutils_build = custom_target('mathutils',
    output:mathutils_out,
    input: mathutils_sources,
    command: [f2py, '-c', '@INPUT@', '--backend=meson', '--quiet', '-m', 'mathutils'],
    env: fc_env,
    install: true,
    install_dir: module_install_dir)

gravity_sources = [
    'tstrippy/src/mathutils.f90',
    'tstrippy/src/sphericalharmonicsbfe.f90',
    'tstrippy/src/besselbfe.f90',
    'tstrippy/src/gravity.f90'
]
gravity_out = 'gravity' + ext_suffix
gravity_build = custom_target('gravity',
    output: gravity_out,
    input: gravity_sources,
    command: [f2py, '-c', '@INPUT@', '--backend=meson', '--quiet', '-m', 'gravity'],
    env: fc_env,
    install: true,
    install_dir: module_install_dir
)

simulator_sources = [
    'tstrippy/src/mathutils.f90',
    'tstrippy/src/sphericalharmonicsbfe.f90',
    'tstrippy/src/besselbfe.f90',
    'tstrippy/src/gravity.f90',
    'tstrippy/src/simulator.f90'
]
simulator_out = 'simulator' + ext_suffix
simulator_build = custom_target('simulator',
    output: simulator_out,
    input: simulator_sources,
    command: [f2py, '-c', '@INPUT@', '--backend=meson', '--quiet', '-m', 'simulator'],
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
        'tstrippy/io/__init__.py',
        'tstrippy/io/mwgcs.py',    
        'tstrippy/io/milkyway_models.py',
        'tstrippy/io/reference_frames.py',
        'tstrippy/io/write_simulation_hdf5.py',
    ),
    subdir: python_install_subdir / 'io'
)

py.install_sources(
    files(
        'tstrippy/code/__init__.py',
        'tstrippy/code/bfe.py',
        'tstrippy/code/sampling.py',
        'tstrippy/code/orbits.py',
    ),
    subdir: python_install_subdir / 'code'
)

# And data files
install_data(
    files('tstrippy/data/mwgcs-kinematics',
          'tstrippy/data/mwgcs-structural',
          'tstrippy/data/MWrefframeFerrone2023.yaml',
          'tstrippy/data/pouliasis2017pii.yaml',
          'tstrippy/data/unit_basis.yaml'),
    install_dir: python_install_subdir / 'data'
)