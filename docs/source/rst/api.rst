Python API Documentation
========================

broken into modules...

Pre-Processor Modules
---------------------

blah blah text

.. _label-sections-module:

*sections* Module
^^^^^^^^^^^^^^^^^

blah blah text

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
""""""""""""""""
..  autoclass:: sectionproperties.pre.sections.Chs
    :show-inheritance:

Rhs Class
""""""""""""""""
..  autoclass:: sectionproperties.pre.sections.Rhs
    :show-inheritance:

ISection Class
""""""""""""""
  ..  autoclass:: sectionproperties.pre.sections.ISection
      :show-inheritance:

PfcSection Class
""""""""""""""""
  ..  autoclass:: sectionproperties.pre.sections.PfcSection
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

blah blah text

..  autofunction:: sectionproperties.pre.pre.create_mesh


Analysis Modules
----------------

blah blah text

*cross_section* Module
^^^^^^^^^^^^^^^^^^^^^^

blah blah text

CrossSection Class
""""""""""""""""""

  ..  autoclass:: sectionproperties.analysis.cross_section.CrossSection
      :show-inheritance:
      :members:

SectionProperties Class
"""""""""""""""""""""""

..  autoclass:: sectionproperties.analysis.cross_section.SectionProperties
    :show-inheritance:
    :members:

StressResult Class
""""""""""""""""""

..  autoclass:: sectionproperties.analysis.cross_section.StressResult
    :show-inheritance:
    :members:

*fea* Module
^^^^^^^^^^^^

blah blah text

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

blah blah text

solver Functions
""""""""""""""""

..  autofunction:: sectionproperties.analysis.solver.solve_cgs
..  autofunction:: sectionproperties.analysis.solver.solve_cgs_lagrange
..  autofunction:: sectionproperties.analysis.solver.solve_direct
..  autofunction:: sectionproperties.analysis.solver.solve_direct_lagrange

Post-Processor Modules
----------------------

blah blah text

*post* Module
^^^^^^^^^^^^^

blah blah text

post Functions
""""""""""""""

..  autofunction:: sectionproperties.post.post.setup_plot
..  autofunction:: sectionproperties.post.post.finish_plot
..  autofunction:: sectionproperties.post.post.print_results
