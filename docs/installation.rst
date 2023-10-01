.. _label-installation:

Installation
============

These instructions will get you a copy of ``sectionproperties`` up and running on your
machine. You will need a working copy of python 3.9, 3.10 or 3.11 to get started.

Installing ``sectionproperties``
--------------------------------

``sectionproperties`` uses `shapely <https://github.com/shapely/shapely>`_ to prepare
the cross-section geometry and `triangle <https://github.com/drufat/triangle>`_ to
efficiently generate a conforming triangular mesh.
`numpy <https://github.com/numpy/numpy>`_ and `scipy <https://github.com/scipy/scipy>`_
are used to aid finite element computations, while
`matplotlib <https://github.com/matplotlib/matplotlib>`_ and
`rich <https://github.com/Textualize/rich>`_ are used for post-processing.
Finally, `click <https://github.com/pallets/click>`_ is used to power the
``sectionproperties`` CLI.

``sectionproperties`` and all of its dependencies can be installed through the python
package index:

.. code-block:: shell

    pip install sectionproperties

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

Note that the ``rhino`` option only supports python ``3.9`` due to incomplete wheel
coverage of ``rhino3dm``.
