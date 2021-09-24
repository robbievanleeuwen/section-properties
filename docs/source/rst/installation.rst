Installation
============

These instructions will get you a copy of *sectionproperties* up and running on
your local machine. You will need a working copy of python 3.8 on your machine.

Installing *sectionproperties*
------------------------------

*sectionproperties* uses *meshpy* to efficiently generate a conforming triangular
mesh in order to perform a finite element analysis of the structural cross-section.
The installation procedure for *meshpy* depends on your local machine.

UNIX (MacOS/Linux)
^^^^^^^^^^^^^^^^^^

*sectionproperties* and all of its dependencies can be installed through the
python package index::

  $ pip install sectionproperties

If you have any issues installing *meshpy*, refer to the installation instructions
on the `MeshPy GitHub page <https://github.com/inducer/meshpy>`_ or the
`MeshPy documentation <https://documen.tician.de/meshpy/installation.html>`_.

If you have any issues installing *shapely*, refer to the installation instructions
on the `Shapely GitHub page <https://github.com/Toblerity/Shapely>`_ or the
`Shapely documentation <https://shapely.readthedocs.io/en/stable/manual.html>`_.

Windows
^^^^^^^

If you use conda, you can install *meshpy* and *shapely* directly::

  $ conda install -c conda-forge meshpy shapely

Alternatively, you can install both *meshpy* and shapely by downloading the appropriate
`installation wheels <https://www.lfd.uci.edu/~gohlke/pythonlibs>`_.

Navigate to the location of the downloaded wheel and install using pip::

  $ cd Downloads
  $ pip install <wheel file names go here>

Once *meshpy* and *shapely* have been installed, the rest of the *sectionproperties* package can
be installed using the python package index::

  $ pip install sectionproperties

Testing the Installation
------------------------

Python *pytest* modules are located in the *sectionproperties.tests* package.
To see if your installation is working correctly, install `pytest` and run the
following test::

  $ pytest --pyargs sectionproperties
