project(
    'tstrippy', 
    ['fortran'], 
    version:'0.0.1',
    license:'MIT',
    default_options:[
        'python.install_env=auto',
    ] )
pkg_name = 'tstrippy'
pymod = import('python')
py = pymod.find_installation(required:true,pure:false)

install_dir = py.get_install_dir()
python_install_subdir = install_dir / pkg_name
message('install dir: ', install_dir)

py.install_sources(
    files(
        'tstrippy/__init__.py',
    ),
    subdir: 'tstrippy' 
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
    subdir:  'tstrippy/code'
)
py.install_sources(
    files(
        'tstrippy/io/__init__.py',
        'tstrippy/io/mwgcs.py',    
        'tstrippy/io/milkyway_models.py',
        'tstrippy/io/reference_frames.py',
        'tstrippy/io/write_simulation_hdf5.py',
    ),
    subdir: 'tstrippy/io'
)
# installing data is native meson
# so we have to specify the full path, there is no "subdir" argument
datadir=python_install_subdir/'data'
install_data(
    files(
        'tstrippy/data/2023-03-28-merged.fits',
        'tstrippy/data/README_baumgardt_catalog.md',
        'tstrippy/data/mwgcs-kinematics',
        'tstrippy/data/mwgcs-structural',
        'tstrippy/data/unit_basis.yaml',
    ),
    install_dir:datadir
)
install_subdir(
    'tstrippy/data/MWmodels',
    install_dir:datadir
)
install_subdir(
    'tstrippy/data/MWreferenceframes',
    install_dir:datadir
)


## START BUILDING FORTRAN CUSOTM TARGETS
# get f2py from environment
ext_suffix_code = '''
import sysconfig
print(sysconfig.get_config_var('EXT_SUFFIX'))
'''
ext_suffix = run_command(py, '-c', ext_suffix_code, check: true).stdout().strip()
f2py_code = '''
import sys, os, numpy
f2py_path = os.path.join(os.path.dirname(numpy.__file__), 'f2py')
print(f2py_path)
'''
f2py_result = run_command(py, '-c', f2py_code, check:true).stdout().strip()
f2py = [py,f2py_result ]
# build mathutils
mathutils_sources =[
    'tstrippy/src/mathutils.f90',
]
mathutils_out = 'mathutils' + ext_suffix
mathutils_build = custom_target('mathutils',
    output:mathutils_out,
    input:mathutils_sources,
    command: [f2py, '-c', '@INPUT@', '--backend=meson', '--quiet', '-m', 'mathutils'],
    install:true,
    install_dir:python_install_subdir)

# NOTE. 
#   Executing python code is necessary to obtain the correct build
#   and install information for running f2py. Meson doesn't have 
#   build in options 
# get the python from the current conda environment
# code_site_packages = '''
# import sysconfig
# print(sysconfig.get_path('purelib'))
# '''
# # e.g. .cpython-312-darwin.so
# code_extension_suffix = '''
# import sysconfig
# print(sysconfig.get_config_var('EXT_SUFFIX'))
# '''
# # just build the python code first 
