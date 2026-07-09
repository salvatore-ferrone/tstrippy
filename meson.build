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
