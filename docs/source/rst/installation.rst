Installation
============

These instructions will get you a copy of *sectionproperties* up and running on
your local machine. You will need a working copy of python 3.8, 3.9 or 3.10 on your
machine.

Installing *sectionproperties*
------------------------------

*sectionproperties* uses `shapely <https://github.com/shapely/shapely>`_ to prepare the
cross-section geometry and `triangle <https://github.com/drufat/triangle>`_ to efficiently
generate a conforming triangular mesh in order to perform a finite element analysis of the
structural cross-section.

*sectionproperties* and all of its dependencies can be installed through the python package index:

.. code-block:: console

  pip install sectionproperties

Note that dependencies required for importing from rhino files are not included by default.
To obtain these dependencies, install using the *rhino* option:

.. code-block:: console

  pip install sectionproperties[rhino]


Testing the Installation
------------------------

Python *pytest* modules are located in the *sectionproperties.tests* package.
To see if your installation is working correctly, install `pytest` and run the
following test:

.. code-block:: console

  pytest --pyargs sectionproperties
