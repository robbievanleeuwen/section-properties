Examples
========

The following examples showcase the features of ``sectionproperties`` in an
instructional manner, while also highlighting several validation tests and advanced
applications. The below examples are available as jupyter notebooks and can be
downloaded
`here <https://github.com/robbievanleeuwen/section-properties/tree/master/docs/examples>`_
(click on the file you would like to download, then click the download icon at the top
of the file).

To run these notebooks you must first install `jupyter notebook <https://jupyter.org/>`_
by running ``pip install notebook`` in the same virtual environment in which
``sectionproperties`` is installed. Next, navigate to the location of the downloaded
examples and run ``jupyter notebook`` to open jupyter in the browser. Finally, double
click an example file to open the notebook, and execute each cell by clicking the play
button. More information on jupyter notebooks can be found
`here <https://docs.jupyter.org/en/latest/>`_. Don't be afraid to
`raise an issue <https://github.com/robbievanleeuwen/section-properties/issues>`_ or
`post in the discussions page <https://github.com/robbievanleeuwen/section-properties/discussions>`_
if you have trouble running any of the examples!

Geometry
--------

.. nbgallery::
    :name: geometry-gallery
    :maxdepth: 1

    examples/geometry/geometry_coordinates
    examples/geometry/section_library
    examples/geometry/geometry_manipulation
    examples/geometry/geometry_cad
    examples/geometry/advanced_geometry
    examples/geometry/create_mesh

Materials
---------

.. nbgallery::
    :name: materials-gallery
    :maxdepth: 1

    examples/materials/assign_materials
    examples/materials/composite_analysis

Analysis
--------

.. nbgallery::
    :name: analysis-gallery
    :maxdepth: 1

    examples/analysis/geometric_analysis
    examples/analysis/plastic_analysis
    examples/analysis/warping_analysis
    examples/analysis/frame_analysis
    examples/analysis/stress_analysis

Results
-------

.. nbgallery::
    :name: results-gallery
    :maxdepth: 1

    examples/results/display_results
    examples/results/get_results
    examples/results/plot_centroids
    examples/results/plot_stress
    examples/results/get_stress
    examples/results/export_fibre_section

Validation
----------

.. nbgallery::
    :name: validation-gallery
    :maxdepth: 1

    examples/validation/pilkey_channel
    examples/validation/pilkey_arc
    examples/validation/pilkey_composite
    examples/validation/peery

Advanced
--------

.. nbgallery::
    :name: advanced-gallery
    :maxdepth: 1

    examples/advanced/advanced_plot
    examples/advanced/rectangle_torsion
    examples/advanced/trapezoidal_torsion
