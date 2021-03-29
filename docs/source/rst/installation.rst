Installation
============

These instructions will get you a copy of *sectionproperties* up and running on
your local machine. You will need a working copy of python>=3.5 on your machine.

Installing *sectionproperties*
------------------------------

*sectionproperties* uses *meshpy* to efficiently generate a conforming triangular
mesh in order to perform a finite element analysis of the structural cross-section.
The installation procedure for of *meshpy* depends on your local machine.

UNIX (MacOS/Linux)
^^^^^^^^^^^^^^^^^^

*sectionproperties* and all of its dependencies can be installed through the
python package index::

  $ pip install sectionproperties

If you have any issues installing *meshpy*, refer to the installation instructions
on its `github page
<https://github.com/inducer/meshpy>`_ or its
`documentation
<https://documen.tician.de/meshpy/installation.html>`_.

Windows
^^^^^^^

Install *meshpy* by downloading the appropriate `installation wheel
<https://www.lfd.uci.edu/~gohlke/pythonlibs/#meshpy>`_.

Navigate to the location of the downloaded wheel and install using pip::

  $ cd Downloads
  $ pip install MeshPy‑2018.2.1‑cp36‑cp36m‑win_amd64.whl

Once *meshpy* has been installed, the rest of the *sectionproperties* package can
be installed using the python package index::

  $ pip install sectionproperties
