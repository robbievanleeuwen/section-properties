.. _label-installation:

Installation
============

These instructions will get you a copy of ``sectionproperties`` up and running on your
machine. You will need a working copy of python 3.10, 3.11 or 3.12 to get started.

Installing ``sectionproperties``
--------------------------------

``sectionproperties`` uses `shapely <https://github.com/shapely/shapely>`_ to prepare
the cross-section geometry and `CyTriangle <https://github.com/m-clare/cytriangle>`_ to
efficiently generate a conforming triangular mesh.
`numpy <https://github.com/numpy/numpy>`_ and `scipy <https://github.com/scipy/scipy>`_
are used to aid finite element computations, while
`matplotlib <https://github.com/matplotlib/matplotlib>`_ and
`rich <https://github.com/Textualize/rich>`_ are used for post-processing.

``sectionproperties`` and all of its dependencies can be installed through the python
package index:

.. code-block:: shell

    pip install sectionproperties

Installing ``Numba``
--------------------

``Numba`` translates a subset of Python and NumPy code into fast machine code, allowing
algorithms to approach the speeds of C. The speed of several ``sectionproperties``
analysis functions have been enhanced with `numba <https://github.com/numba/numba>`_.
To take advantage of this increase in performance you can install ``numba`` alongside
``sectionproperties`` with:

.. code-block:: shell

    pip install sectionproperties[numba]

Installing ``PARDISO`` Solver
-----------------------------

The default sparse solver used in ``scipy`` is ``SuperLU``.
It performs okay for small matrices but appears to be slower for larger matrices. The
``PARDISO`` solver is a much faster alternative
(see `pypardiso <https://github.com/haasad/PyPardisoProject>`_), but it requires the
installation of the ``MKL`` library, which takes a lot of disk space. Note that this
library is only available for Linux and Windows systems.

If you do not have a disk space constraint, you can install the ``PARDISO`` solver with:

.. code-block:: shell

    pip install sectionproperties[pardiso]

Installing CAD Modules
----------------------

The dependencies used to import from ``.dxf`` and ``.3dm`` (rhino) files are not
included by default in the base installation.
`cad-to-shapely <https://github.com/aegis1980/cad-to-shapely>`_ is used to import
``.dxf`` files, while
`rhino-shapely-interop <https://github.com/normanrichardson/rhino_shapely_interop>`_ is
used to import ``.3dm`` files.

To install ``sectionproperties`` with the above functionality, use the ``dxf`` and/or
``rhino`` options:

.. code-block:: shell

    pip install sectionproperties[dxf]
    pip install sectionproperties[rhino]
