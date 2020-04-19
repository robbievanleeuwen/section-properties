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

PolygonSection Class
""""""""""""""""""""
  ..  autoclass:: sectionproperties.pre.sections.PolygonSection
      :show-inheritance:

BoxGirderSection Class
""""""""""""""""""""""
  ..  autoclass:: sectionproperties.pre.sections.BoxGirderSection
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


*offset* Module
^^^^^^^^^^^^^^^

..  autofunction:: sectionproperties.pre.offset.offset_perimeter


*nastran_sections* Module
^^^^^^^^^^^^^^^^^^^^^^^^^
This module contains cross-sections as defined by Nastran and Nastran-based programs,
such as MYSTRAN and ASTROS.

BARSection Class
""""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.BARSection
    :show-inheritance:
    :members:

BOXSection Class
""""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.BOXSection
    :show-inheritance:
    :members:

BOX1Section Class
"""""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.BOX1Section
    :show-inheritance:
    :members:

CHANSection Class
"""""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.CHANSection
    :show-inheritance:
    :members:

CHAN1Section Class
""""""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.CHAN1Section
    :show-inheritance:
    :members:

CHAN2Section Class
""""""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.CHAN2Section
    :show-inheritance:
    :members:

CROSSSection Class
"""""""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.CROSSSection
    :show-inheritance:
    :members:

DBOXSection Class
"""""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.DBOXSection
    :show-inheritance:
    :members:

FCROSSSection Class
"""""""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.FCROSSSection
    :show-inheritance:
    :members:

GBOXSection Class
"""""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.GBOXSection
    :show-inheritance:
    :members:

HSection Class
""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.HSection
    :show-inheritance:
    :members:

HATSection Class
""""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.HATSection
    :show-inheritance:
    :members:

HAT1Section Class
"""""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.HAT1Section
    :show-inheritance:
    :members:

HEXASection Class
"""""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.HEXASection
    :show-inheritance:
    :members:

NISection Class
"""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.NISection
    :show-inheritance:
    :members:

I1Section Class
"""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.I1Section
    :show-inheritance:
    :members:

LSection Class
""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.LSection
    :show-inheritance:
    :members:

RODSection Class
""""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.RODSection
    :show-inheritance:
    :members:

TSection Class
""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.TSection
    :show-inheritance:
    :members:

T1Section Class
"""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.T1Section
    :show-inheritance:
    :members:

T2Section Class
"""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.T2Section
    :show-inheritance:
    :members:

TUBESection Class
"""""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.TUBESection
    :show-inheritance:
    :members:

TUBE2Section Class
""""""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.TUBE2Section
    :show-inheritance:
    :members:

ZSection Class
""""""""""""""
..  autoclass:: sectionproperties.pre.nastran_sections.ZSection
    :show-inheritance:
    :members:

References
""""""""""
.. [1]  MSC Nastran Quick Reference Guide 2012,
        PBEAML - Simple Beam Cross-Section Property, pp. 2890-2894
        https://simcompanion.mscsoftware.com/infocenter/index?page=content&id=DOC10351
.. [2]  Siemens NX Nastran 12 Quick Reference Guide,
        PBEAML, pp. 16-59 - 16-62
        https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/12/help/tdoc/en_US/pdf/QRG.pdf
.. [3]  AutoDesk Nastran Online Documentation, Nastran Reference Guide,
        Section 4 - Bulk Data, PBEAML
        http://help.autodesk.com/view/NSTRN/2018/ENU/?guid=GUID-B7044BA7-3C26-49DA-9EE7-DA7505FD4B2C
.. [4]  Users Reference Manual for the MYSTRAN General Purpose Finite Element Structural Analysis Computer Program,
        Jan. 2019, Section 6.4.1.53 - PBARL, pp. 154-156
        https://www.mystran.com/Executable/MYSTRAN-Users-Manual.pdf
.. [5]  Astros Enhancements - Volume III - Astros Theoretical Manual,
        Section 5.1.3.2, pp. 56
        https://apps.dtic.mil/dtic/tr/fulltext/u2/a308134.pdf

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
