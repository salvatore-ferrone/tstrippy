from numpy.distutils.core import setup, Extension
import setuptools
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

# my fortran module


constands_module = Extension('tstrippy.lib.constants',
                              sources=[
                                  'tstrippy/src/constants.f90',
                                  ])

potentials_module = Extension('tstrippy.lib.potentials',
                              sources=[
                                  'tstrippy/src/constants.f90',
                                  'tstrippy/src/potentials.f90',
                                  ])

integrator_module = Extension('tstrippy.lib.integrator',
                              sources=[
                                  'tstrippy/src/constants.f90',
                                  'tstrippy/src/potentials.f90',
                                  'tstrippy/src/galacticbar.f90',
                                  'tstrippy/src/perturbers.f90',
                                  'tstrippy/src/hostperturber.f90',
                                  'tstrippy/src/integrator.f90', 
                                  ])



setup(
    name="tstrippy",
    version="0.0.1",
    author="Salvatore Ferrone",
    author_email="salvatore.ferrone.1996@gmail.com",
    description="A small example package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/salvatore-ferrone",
    packages=setuptools.find_packages(include=['tstrippy', 'tstrippy.*']),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Mac only (for now)",
    ],
    python_requires='>=3.6',
    install_requires=['numpy',
                      "h5py",
                      'astropy',],
    ext_modules=[constands_module,potentials_module,integrator_module]
)