Welcome to tstrippy's documentation
===================================

tstrippy is a research code for modeling stellar stream formation in the Milky Way,
with an initial focus on stars escaping from globular clusters. It provides
gravitational potential models and numerical tools for integrating Hamiltonian
equations of motion.

This package was developed by `Salvatore Ferrone <https://salvatore-ferrone.github.io/>`_
for the simulations presented in the
`PhD thesis <https://salvatore-ferrone.github.io/phd-thesis-ferrone/>`_.
It is not intended to replace mature community packages such as
`galpy <https://docs.galpy.org/en/v1.11.2/>`_,
`gala <https://gala.adrian.pw/en/latest/#>`_, or
`Agama <https://agama.software/>`_. Instead, it is a focused tool built for a
specific scientific problem.

If you use this package and have questions, feel free to reach out:
salvatore.ferrone.199@gmail.com.
Precomputed predictions for escaped tidal debris from Milky Way globular clusters
are available `here <https://etidal-project.obspm.fr/simulations.html>`_.

Documentation roadmap
---------------------

This documentation is organized into three sections:

1. Package basics:
   loading included data products (currently the globular cluster catalog) and
   working with the built-in gravitational potential models.
2. Orbit integration in the Galaxy:
   solving for trajectories of point masses in different Galactic potentials.
3. Stream disruption workflows:
   simulating stream formation/disruption and validating that the method behaves
   as expected.

.. toctree::
   :maxdepth: 1
   :caption: Installation

   installation.rst

.. toctree::
   :maxdepth: 1
   :caption: Package basics

   baumgardt_catalog.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Orbit integration in the Milky Way

   galactic_orbits.ipynb
   reverse_integrability_vanilla.ipynb
   bar.ipynb
   reverse_integrability_bar.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Stream disruption examples

   bar_stream_example.ipynb
   reverse_integrability_stream.ipynb