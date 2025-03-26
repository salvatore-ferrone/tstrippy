project('hello_project', 'fortran', version: '0.1', default_options: ['buildtype=release'])

# Import Python module
py3 = import('python').find_installation('python3', required: true)

# Define the Python extension module
hello_module = py3.extension_module(
  'hello',  # Name of the Python module
  sources: 'development/hello.f90',
  install: true
)