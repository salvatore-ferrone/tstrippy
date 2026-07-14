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

# the snippit below is preferd to: install_dir = py.get_install_dir()
install_dir_code='''import sysconfig; print(sysconfig.get_path('platlib'));'''
install_dir = run_command(py, '-c', install_dir_code, check:true).stdout().strip()
python_install_subdir = install_dir / pkg_name

message('install dir: ', install_dir)

# install the python part of the project
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
# make sure we're using f2py in the conda environment and not the system-wide version
f2py_code = '''
import sys, os, numpy
f2py_path = os.path.join(os.path.dirname(numpy.__file__), 'f2py')
print(f2py_path)
'''
f2py_result = run_command(py, '-c', f2py_code, check:true).stdout().strip()
f2py = [py,f2py_result] # essenctiall: python -m numpy.f2py

# 1. Add Agama
# 2. implement the option for Agama
# 3. implement the stub file

agama_so_fname = 'agama.so'
agama_so_path = install_dir / 'agama' / agama_so_fname


# Environment for targets that include agamabackend.f90.
# -undefined dynamic_lookup (macOS) defers resolution of Agama symbols to
# runtime so the .so loads cleanly even when Agama is not linked in.



# fc_env = environment()
# fc_env.set('FC', 'gfortran')
# fc_env.set('FFLAGS', '-Wno-unused-variable -fPIC')
# fc_env.set('CFLAGS', '-Wno-misleading-indentation')  # Suppress misleading indentation warning
# fc_env.set('LDFLAGS', agama_link_flags.strip())

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
    command: [
        f2py,
        '-c',
        '@INPUT@',
        '--backend=meson',
        '--quiet',
        '-m',
        'gravity',
        agama_so_path,
    ],
    # env: fc_env,
    install:true,
    install_dir:python_install_subdir
)


simulator_sources = [
    'tstrippy/src/mathutils.f90',
    'tstrippy/src/sphericalharmonicsbfe.f90',
    'tstrippy/src/besselbfe.f90',
    'tstrippy/src/agamabackend.f90',
    'tstrippy/src/gravity.f90',
    'tstrippy/src/hostcluster.f90',
    'tstrippy/src/simulator.f90',
]
simulator_out = 'simulator' + ext_suffix
simulator_build = custom_target('simulator',
    output: simulator_out,
    input: simulator_sources,
    command: [
        f2py, 
        '-c', 
        '@INPUT@', 
        '--backend=meson', 
        '--quiet', 
        '-m', 
        'simulator',
        agama_so_path],
    # env: fc_env,
    install: true,
    install_dir: python_install_subdir
)


# now update the dynamic link for macOS. 
# this is required because at run time, simulator and gravity will look for agama.so
#   only in /usr/lib/. We have to add the symbolic link to where agama.so is located

if host_machine.system()=='darwin'
    install_name_tool = find_program('install_name_tool', required:true)
    meson.add_install_script(
        install_name_tool,
        '-change',
        '@executable_path/'+agama_so_fname,
        agama_so_path,
        install_dir / pkg_name / simulator_out
    )
    meson.add_install_script(
        install_name_tool,
        '-change',
        '@executable_path/'+agama_so_fname,
        agama_so_path,
        install_dir / pkg_name / gravity_out
    )
endif