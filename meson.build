project('tstrippy', ['c'],
    version: '0.0.1',
    license: 'MIT'
)

# Define package name variable to avoid repetition
pkg_name = 'tstrippy'

# Find Python and f2py
py = import('python').find_installation(required:true, pure:false)

# Get the Python site-packages directory
code_site_packages = '''
import sysconfig
print(sysconfig.get_path('purelib'))
'''
# Get Python extension suffix for module naming
ext_suffix_code = '''
import sysconfig
print(sysconfig.get_config_var('EXT_SUFFIX'))
'''
conda_prefix_code = '''
import sys
print(sys.prefix)
'''

install_dir = run_command(py, '-c', code_site_packages, check: true).stdout().strip()
ext_suffix = run_command(py, '-c', ext_suffix_code, check: true).stdout().strip()
python_prefix = run_command(py, '-c', python_prefix_code, check: true).stdout().strip()
python_lib_dir = python_prefix / 'lib'
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
fc_args = ['-fPIC']

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
fs = import('fs')

# Optional Agama wiring
agama_opt = get_option('agama')
agama_dir = meson.project_source_root() / 'external' / 'agama'
agama_a = agama_dir / 'agama.a'
agama_so = agama_dir / 'agama.so'
openblas_so = python_lib_dir / 'libopenblas.0.dylib'
gsl_a = python_lib_dir / 'libgsl.a'
gslcblas_a = python_lib_dir / 'libgslcblas.a'
libomp_so = python_lib_dir / 'libomp.dylib'
agama_enabled = false
agama_link_flags = ''

if not agama_opt.disabled()
    if fs.exists(agama_a)
        agama_enabled = true
        if fs.exists(openblas_so)
            agama_link_flags = '-Wl,-rpath,' + python_lib_dir + ' ' + agama_a + ' ' + gsl_a + ' ' + gslcblas_a + ' ' + openblas_so + ' ' + libomp_so + ' -lc++'
        else
            agama_link_flags = '-Wl,-rpath,' + python_lib_dir + ' ' + agama_a + ' ' + gsl_a + ' ' + gslcblas_a + ' -lopenblas ' + libomp_so + ' -lc++'
        endif
        message('Agama backend enabled via ' + agama_a)
    elif fs.exists(agama_so)
        agama_enabled = true
        agama_link_flags = '-Wl,-rpath,' + agama_dir + ' ' + agama_so
        message('Agama backend enabled via ' + agama_so)
    elif agama_opt.enabled()
        error('Agama backend explicitly enabled, but external/agama/agama.a or external/agama/agama.so was not found. Build Agama first with make lib in external/agama.')
    else
        warning('Agama backend not found at external/agama/agama.so; continuing without Agama support')
    endif
endif


# Setup environment for f2py
base_ldflags = ''

# Suppress compact unwind warnings from macOS linker
if host_machine.system() == 'darwin'
    # For macOS, add linker flags to suppress unwind warnings
    base_ldflags = base_ldflags + ' -Wl,-no_compact_unwind'
endif

fc_env_base = environment()
fc_env_base.set('FC', 'gfortran')
fc_env_base.set('FFLAGS', '-Wno-unused-variable -fPIC')
fc_env_base.set('CFLAGS', '-Wno-misleading-indentation')  # Suppress misleading indentation warning
if base_ldflags.strip() != ''
    fc_env_base.set('LDFLAGS', base_ldflags.strip())
endif

fc_env_agama = environment()
fc_env_agama.set('FC', 'gfortran')
fc_env_agama.set('FFLAGS', '-Wno-unused-variable -fPIC')
fc_env_agama.set('CFLAGS', '-Wno-misleading-indentation')  # Suppress misleading indentation warning

agama_ldflags = base_ldflags
if agama_enabled
    agama_ldflags = agama_ldflags + ' ' + agama_link_flags
endif
if agama_ldflags.strip() != ''
    fc_env_agama.set('LDFLAGS', agama_ldflags.strip())
endif

agamabackend_sources = []
if agama_enabled
    agamabackend_sources += ['tstrippy/src/agamabackend.f90']
endif

mathutils_sources = [
    'tstrippy/src/mathutils.f90',
]
mathutils_out = 'mathutils' + ext_suffix
mathutils_build = custom_target('mathutils',
    output:mathutils_out,
    input: mathutils_sources,
    command: [f2py, '-c', '@INPUT@', '--backend=meson', '--quiet', '-m', 'mathutils'],
    env: fc_env_base,
    install: true,
    install_dir: module_install_dir)

gravity_sources = [
    'tstrippy/src/mathutils.f90',
    'tstrippy/src/sphericalharmonicsbfe.f90',
    'tstrippy/src/besselbfe.f90',
    'tstrippy/src/agamabackend.f90',
    'tstrippy/src/gravity.f90',
]
gravity_out = 'gravity' + ext_suffix
gravity_build = custom_target('gravity',
    output: gravity_out,
    input: gravity_sources,
    command: [f2py, '-c', '@INPUT@', '--backend=meson', '--quiet', '-m', 'gravity'],
    env: fc_env_agama,
    install: true,
    install_dir: module_install_dir
)

simulator_sources = [
    'tstrippy/src/mathutils.f90',
    'tstrippy/src/sphericalharmonicsbfe.f90',
    'tstrippy/src/besselbfe.f90',
    'tstrippy/src/agamabackend.f90',
    'tstrippy/src/gravity.f90',
    'tstrippy/src/hostcluster.f90',
    'tstrippy/src/simulator.f90'
]
simulator_out = 'simulator' + ext_suffix
simulator_build = custom_target('simulator',
    output: simulator_out,
    input: simulator_sources,
    command: [f2py, '-c', '@INPUT@', '--backend=meson', '--quiet', '-m', 'simulator'],
    env: fc_env_agama,
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
        'tstrippy/code/analysis.py',
        'tstrippy/code/workflows.py',
    ),
    subdir: python_install_subdir / 'code'
)

# for listing files 

# colelct all files in the MW referenceframes 


# And data files
# install_data(
#     files('tstrippy/data/mwgcs-kinematics',
#           'tstrippy/data/mwgcs-structural',
#           'tstrippy/data/MWreferenceframes/MWrefframeFerrone2023.yaml',
#           'tstrippy/data/MWmodels/pouliasis2017pii.yaml',
#           'tstrippy/data/MWmodels/ibata2024mostprob.yaml',
#           'tstrippy/data/unit_basis.yaml'),
#     install_dir: python_install_subdir / 'data'
# )

# And data files
# Keep top-level data files explicit, and install model/reference-frame
# directories as subtrees so new files are picked up automatically.
install_data(
    files(
        'tstrippy/data/2023-03-28-merged.fits',
        'tstrippy/data/README_baumgardt_catalog.md',
        'tstrippy/data/mwgcs-kinematics',
        'tstrippy/data/mwgcs-structural',
        'tstrippy/data/unit_basis.yaml',
    ),
    install_dir: python_install_subdir / 'data'
)

install_subdir(
    'tstrippy/data/MWmodels',
    install_dir: python_install_subdir / 'data'
)

install_subdir(
    'tstrippy/data/MWreferenceframes',
    install_dir: python_install_subdir / 'data'
)