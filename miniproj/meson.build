project('miniproj', ['fortran'], version:'0.0.1',default_options:['python.install_env=auto'])
pkg_name = 'miniproj'
pymod = import('python')
py = pymod.find_installation(required:true,pure:false)
python_install_dir = py.get_install_dir()
python_install_subdir = python_install_dir / pkg_name
message('install dir: ', python_install_dir)

py.install_sources(
    files(
        join_paths(pkg_name,'__init__.py')
    ),
    subdir: pkg_name 
)

##############################################################################################
######################## PREPARE ARGUMENTS FOR FORTRAN CUSTOM TARGETS ########################
##############################################################################################
ext_suffix_code = '''
import sysconfig
print(sysconfig.get_config_var('EXT_SUFFIX'))
'''
ext_suffix = run_command(py, '-c', ext_suffix_code, check: true).stdout().strip()
# on macOS: EXT_SUFFIX=.cpython-312-darwin.so
# on linux: EXT_SUFFIX=.cpython-312-x86_64-linux-gnu.so
#   if python3.12, for example. 

f2py_code = '''
import sys, os, numpy
f2py_path = os.path.join(os.path.dirname(numpy.__file__), 'f2py')
print(f2py_path)
'''
f2py_result = run_command(py, '-c', f2py_code, check:true).stdout().strip()
f2py = [py,f2py_result] # essentially: python -m numpy.f2py

fc_env = environment()
fc_env.set('FC', 'gfortran')
fc_env.set('FFLAGS', '-Wno-unused-variable -fPIC')
fc_env.set('CFLAGS', '-Wno-misleading-indentation -Wno-unused-but-set-variable')  


# prepare the arguments for building the fotran shared object files 
simulator_out = 'simulator' + ext_suffix
simulator_sources =[pkg_name + '/src/simulator.f90', pkg_name + '/src/particles.f90',]
simulator_command = [f2py, '-c', '@INPUT@', '--backend=meson', '--quiet', '-m', 'simulator']

simbuild = custom_target('simulator',
    output:simulator_out,
    input:simulator_sources,
    command:simulator_command,
    env:fc_env,
    install:true,
    install_dir:python_install_subdir)