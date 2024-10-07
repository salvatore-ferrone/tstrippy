Installation instructions
=========================
Dependencies
-------------
``tstrippy`` requires the use of `numpy <https://numpy.org/>`__ version 1.23 due to `f2py` compatibility. There other dependencies as well. To ensure proper function, install the dependencies using the provided `environment.yml` file. This requires conda to be installed. To create a new conda environment with the necessary dependencies, run:


.. code-block:: bash

    conda env create -f environment.yml

This will create a new conda environment called `tstrippy` with all the necessary dependencies. To activate the environment, run:

.. code-block:: bash

    conda activate tstrippy



Installation
-------------
tstrippy cannot be installed via pip. To install tstrippy,  clone or download the source code from `GitHub`_. Navigate to `tstrippy/tstrippy` directory and run the following command:

.. code-block:: bash

    python setup.py install

.. _GitHub: https://github.com/salvatore-ferrone/tstrippy