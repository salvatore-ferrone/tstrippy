from numpy.distutils.core import setup, Extension
import setuptools
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

# my fortran module
ext = Extension('test',sources=['test.f90'])

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
    ext_modules=[ext]
)