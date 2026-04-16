.. image:: _static/logo-light-mode.png
    :class: only-light

.. image:: _static/logo-dark-mode.png
    :class: only-dark

.. toctree::
    :hidden:

    installation
    user_guide
    examples
    api

.. toctree::
    :caption: Development
    :hidden:

    contributing
    Code of Conduct <codeofconduct>
    License <license>
    Changelog <https://github.com/robbievanleeuwen/section-properties/releases>

Documentation
=============

``sectionproperties`` is a python package for the analysis of arbitrary
cross-sections using the finite element method. ``sectionproperties``
can be used to determine section properties to be used in structural
design and visualise cross-sectional stresses resulting from
combinations of applied forces and bending moments.

`Subscribe <http://eepurl.com/dMMUeg>`_ to the ``sectionproperties`` mailing list!

Installation
------------

You can install ``sectionproperties`` via `pip <https://pip.pypa.io/>`_ from
`PyPI <https://pypi.org/>`_:

.. code:: shell

   pip install sectionproperties

See :ref:`label-installation` for more information.

Features
--------

See the complete list of ``sectionproperties`` features :ref:`here<label-features>`.

Quick Start
-----------

Analyse a rectangular cross-section and retrieve some key properties:

.. code:: python

   from sectionproperties.pre.library import rectangular_section
   from sectionproperties.analysis import Section

   # create a 50 x 100 rectangle and mesh it
   geom = rectangular_section(d=100, b=50)
   geom.create_mesh(mesh_sizes=[5])

   # run a geometric analysis
   sec = Section(geometry=geom)
   sec.calculate_geometric_properties()

   # get some results
   area = sec.get_area()
   ixx_c, iyy_c, ixy_c = sec.get_ic()
   print(f"Area = {area:.0f} mm²")
   print(f"Ixx = {ixx_c:.0f} mm⁴, Iyy = {iyy_c:.0f} mm⁴")

.. code:: text

   Area = 5000 mm²
   Ixx = 4166667 mm⁴, Iyy = 1041667 mm⁴

See the :doc:`examples <examples>` for more detailed walkthroughs including composite
sections, warping analysis, and stress calculations.

Contributing
------------

Contributions are very welcome. To learn more, see the
:ref:`Contributor Guide<label-contributing>`.

License
-------

Distributed under the terms of the :doc:`MIT License <license>` ``sectionproperties``
is free and open source software.

Support
-------

Found a bug 🐛, or have a feature request ✨, raise an issue on the
GitHub `issue
tracker <https://github.com/robbievanleeuwen/section-properties/issues>`_.
Alternatively you can get support on the
`discussions <https://github.com/robbievanleeuwen/section-properties/discussions>`_
page.

Disclaimer
----------

``sectionproperties`` is an open source engineering tool that continues to benefit from
the collaboration of many contributors. Although efforts have been made to ensure the
that relevant engineering theories have been correctly implemented, it remains the
user's responsibility to confirm and accept the output. Refer to the
:doc:`License <license>` for clarification of the conditions of use.
