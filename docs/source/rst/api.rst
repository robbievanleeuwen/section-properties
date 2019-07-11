Python API Documentation
========================

Pre-Processor Package
---------------------

.. _label-sections-module:

*sections* Module
^^^^^^^^^^^^^^^^^

Geometry Class
""""""""""""""
..  autoclass:: sectionproperties.pre.sections.Geometry
    :members:

CustomSection Class
"""""""""""""""""""
..  autoclass:: sectionproperties.pre.sections.CustomSection
    :show-inheritance:

RectangularSection Class
""""""""""""""""""""""""
..  autoclass:: sectionproperties.pre.sections.RectangularSection
    :show-inheritance:

CircularSection Class
"""""""""""""""""""""
..  autoclass:: sectionproperties.pre.sections.CircularSection
    :show-inheritance:

Chs Class
"""""""""
..  autoclass:: sectionproperties.pre.sections.Chs
    :show-inheritance:

EllipticalSection Class
"""""""""""""""""""""""
..  autoclass:: sectionproperties.pre.sections.EllipticalSection
    :show-inheritance:

Ehs Class
"""""""""
..  autoclass:: sectionproperties.pre.sections.Ehs
    :show-inheritance:

Rhs Class
"""""""""
..  autoclass:: sectionproperties.pre.sections.Rhs
    :show-inheritance:

ISection Class
""""""""""""""
  ..  autoclass:: sectionproperties.pre.sections.ISection
      :show-inheritance:

MonoISection Class
""""""""""""""""""
  ..  autoclass:: sectionproperties.pre.sections.MonoISection
      :show-inheritance:

TaperedFlangeISection Class
"""""""""""""""""""""""""""
  ..  autoclass:: sectionproperties.pre.sections.TaperedFlangeISection
      :show-inheritance:

PfcSection Class
""""""""""""""""
  ..  autoclass:: sectionproperties.pre.sections.PfcSection
      :show-inheritance:

TaperedFlangeChannel Class
""""""""""""""""""""""""""
  ..  autoclass:: sectionproperties.pre.sections.TaperedFlangeChannel
      :show-inheritance:

TeeSection Class
""""""""""""""""
  ..  autoclass:: sectionproperties.pre.sections.TeeSection
      :show-inheritance:

AngleSection Class
""""""""""""""""""
  ..  autoclass:: sectionproperties.pre.sections.AngleSection
      :show-inheritance:

CeeSection Class
""""""""""""""""
  ..  autoclass:: sectionproperties.pre.sections.CeeSection
      :show-inheritance:

ZedSection Class
""""""""""""""""
  ..  autoclass:: sectionproperties.pre.sections.ZedSection
      :show-inheritance:

CruciformSection Class
""""""""""""""""""""""
  ..  autoclass:: sectionproperties.pre.sections.CruciformSection
      :show-inheritance:

MergedSection Class
"""""""""""""""""""
  ..  autoclass:: sectionproperties.pre.sections.MergedSection
      :show-inheritance:


*pre* Module
^^^^^^^^^^^^

Material Class
""""""""""""""

  ..  autoclass:: sectionproperties.pre.pre.Material
      :show-inheritance:
      :members:

GeometryCleaner Class
"""""""""""""""""""""

  ..  autoclass:: sectionproperties.pre.pre.GeometryCleaner
      :show-inheritance:
      :members:

pre Functions
"""""""""""""

..  autofunction:: sectionproperties.pre.pre.create_mesh


Analysis Package
----------------

*cross_section* Module
^^^^^^^^^^^^^^^^^^^^^^

CrossSection Class
""""""""""""""""""

  ..  autoclass:: sectionproperties.analysis.cross_section.CrossSection
      :show-inheritance:
      :members:

PlasticSection Class
""""""""""""""""""""

  ..  autoclass:: sectionproperties.analysis.cross_section.PlasticSection
      :show-inheritance:
      :members:

StressPost Class
""""""""""""""""

..  autoclass:: sectionproperties.analysis.cross_section.StressPost
    :show-inheritance:
    :members:

MaterialGroup Class
"""""""""""""""""""

  ..  autoclass:: sectionproperties.analysis.cross_section.MaterialGroup
      :show-inheritance:
      :members:

StressResult Class
""""""""""""""""""

..  autoclass:: sectionproperties.analysis.cross_section.StressResult
    :show-inheritance:
    :members:

SectionProperties Class
"""""""""""""""""""""""

..  autoclass:: sectionproperties.analysis.cross_section.SectionProperties
    :show-inheritance:
    :members:

*fea* Module
^^^^^^^^^^^^

Tri6 Class
""""""""""

..  autoclass:: sectionproperties.analysis.fea.Tri6
    :show-inheritance:
    :members:

fea Functions
"""""""""""""

..  autofunction:: sectionproperties.analysis.fea.gauss_points
..  autofunction:: sectionproperties.analysis.fea.shape_function
..  autofunction:: sectionproperties.analysis.fea.extrapolate_to_nodes
..  autofunction:: sectionproperties.analysis.fea.principal_coordinate
..  autofunction:: sectionproperties.analysis.fea.global_coordinate
..  autofunction:: sectionproperties.analysis.fea.point_above_line

*solver* Module
^^^^^^^^^^^^^^^

solver Functions
""""""""""""""""

..  autofunction:: sectionproperties.analysis.solver.solve_cgs
..  autofunction:: sectionproperties.analysis.solver.solve_cgs_lagrange
..  autofunction:: sectionproperties.analysis.solver.solve_direct
..  autofunction:: sectionproperties.analysis.solver.solve_direct_lagrange
..  autofunction:: sectionproperties.analysis.solver.function_timer

Post-Processor Package
----------------------

*post* Module
^^^^^^^^^^^^^

post Functions
""""""""""""""

..  autofunction:: sectionproperties.post.post.setup_plot
..  autofunction:: sectionproperties.post.post.finish_plot
..  autofunction:: sectionproperties.post.post.draw_principal_axis
..  autofunction:: sectionproperties.post.post.print_results
