Installation
============

These instructions will get you a copy of *sectionproperties* up and running on
your local machine. You will need a working copy of python 3.7, 3.8 or 3.9 on your machine.

Installing *sectionproperties*
------------------------------

*sectionproperties* uses `shapely <https://github.com/shapely/shapely>`_ to prepare the
cross-section geometry and `triangle <https://github.com/drufat/triangle>`_ to efficiently
generate a conforming triangular mesh in order to perform a finite element analysis of the
structural cross-section.

*sectionproperties* and all of its dependencies can be installed through the python package index::

  $ pip install sectionproperties

Testing the Installation
------------------------

Python *pytest* modules are located in the *sectionproperties.tests* package.
To see if your installation is working correctly, install `pytest` and run the
following test::

  $ pytest --pyargs sectionproperties
