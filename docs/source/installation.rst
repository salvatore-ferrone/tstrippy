Installation
============

This page describes how to install ``tstrippy`` from source.

Prerequisites
-------------

Before building the package, make sure the following tools are available on your system:

* ``conda`` or ``mamba``
* Python 3.9 to 3.11
* A compatible Fortran compiler such as ``gfortran``

The repository, source code, and issue tracker are available on the
`GitHub repository <https://github.com/salvatore-ferrone/tstrippy>`_.

Clone the repository
--------------------

.. code-block:: bash

    git clone https://github.com/salvatore-ferrone/tstrippy.git
    cd tstrippy

Create the environment
----------------------

The recommended installation path uses the provided conda environment,
which installs the Python dependencies together with the build tools.

.. code-block:: bash

    conda env create -f environment.yml
    conda activate tstrippy

Build the package
-----------------

``tstrippy`` uses ``meson`` as its build system.

.. code-block:: bash

    meson setup builddir
    meson compile -C builddir
    meson install -C builddir

If you prefer, you can use the helper script in the repository root,
which runs the same build steps:

.. code-block:: bash

    ./build.sh

Verify the installation
-----------------------

After the build completes, verify that Python can import the package:

.. code-block:: bash

    python -c "import tstrippy; print(tstrippy.__file__)"

