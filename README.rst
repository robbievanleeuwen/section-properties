sectionproperties
=================

|Build Status| |Documentation Status|

a python package for the analysis of arbitrary cross-sections using the
finite element method written by Robbie van Leeuwen. *sectionproperties*
can be used to determine section properties to be used in structural
design and visualise cross-sectional stresses resulting from
combinations of applied forces and bending moments.

.. raw:: html

    <embed>
      <!-- Begin Mailchimp Signup Form -->
      <link href="//cdn-images.mailchimp.com/embedcode/slim-10_7.css" rel="stylesheet" type="text/css">
      <style type="text/css">
      #mc_embed_signup{background:#fff; clear:left; font:14px Helvetica,Arial,sans-serif; }
      /* Add your own Mailchimp form style overrides in your site stylesheet or in this style block.
       We recommend moving this block and the preceding CSS link to the HEAD of your HTML file. */
      </style>
      <div id="mc_embed_signup">
      <form action="https://github.us19.list-manage.com/subscribe/post?u=541c65ecb1b23522bcf1300db&amp;id=b7a47b4e83" method="post" id="mc-embedded-subscribe-form" name="mc-embedded-subscribe-form" class="validate" target="_blank" novalidate>
      <div id="mc_embed_signup_scroll">
      <label for="mce-EMAIL">Subscribe to the sectionproperties mailing list!</label>
      <input type="email" value="" name="EMAIL" class="email" id="mce-EMAIL" placeholder="email address" required>
      <!-- real people should not fill this in and expect good things - do not remove this or risk form bot signups-->
      <div style="position: absolute; left: -5000px;" aria-hidden="true"><input type="text" name="b_541c65ecb1b23522bcf1300db_b7a47b4e83" tabindex="-1" value=""></div>
      <div class="clear"><input type="submit" value="Subscribe" name="subscribe" id="mc-embedded-subscribe" class="button"></div>
      </div>
      </form>
      </div>
      <!--End mc_embed_signup-->
    </embed>

Installation:
-------------

For more detailed installation instructions, refer to the `documentation`_.

UNIX (MacOS/Linux):
~~~~~~~~~~~~~~~~~~~

::

   $ pip install sectionproperties

Windows
~~~~~~~

Install *meshpy* by downloading the appropriate `installation wheel`_.

Navigate to the location of the downloaded wheel and install using pip:

::

   $ cd Downloads
   $ pip install MeshPy‑2018.2.1‑cp36‑cp36m‑win_amd64.whl

Once *meshpy* has been installed, *sectionproperties* can be installed:

::

   $ pip install sectionproperties

Documentation:
--------------

*sectionproperties* has a fully documented python API which you can find
at https://sectionproperties.readthedocs.io/. To read more about the
theory behind the program, its implementation and some more examples,
check out my blog at https://robbievanleeuwen.github.io/.

Current Capabilities:
---------------------

Pre-Processor:
~~~~~~~~~~~~~~

-  [x] Python API
-  [x] Custom section geometry input
-  [x] Common section geometry generators
-  [x] Multiple geometry merging
-  [x] Geometry cleaning
-  [ ] JSON input file
-  [ ] .dxf import
-  [x] Quadratic triangular mesh generation
-  [x] Composite material properties

Cross-Section Analysis:
~~~~~~~~~~~~~~~~~~~~~~~

-  [x] Global axis geometric section properties:

   -  [x] Area
   -  [x] First moments of area
   -  [x] Second moments of area
   -  [x] Elastic centroid
-  [x] Centroidal axis geometric section properties:
   -  [x] Second moments of area
   -  [x] Elastic section moduli
   -  [ ] Yield moment
   -  [x] Radii of gyration
   -  [x] Plastic centroid
   -  [x] Plastic section moduli
   -  [x] Shape factors
-  [x] Principal axis geometric section properties:
   -  [x] Second moments of area
   -  [x] Elastic section moduli
   -  [ ] Yield moment
   -  [x] Radii of gyration
   -  [x] Plastic centroid
   -  [x] Plastic section moduli
   -  [x] Shape factors
-  [x] Warping section properties:
   -  [x] Torsion constant
   -  [x] Warping constant
-  [x] Shear section properties:
   -  [x] Shear centre (elastic method)
   -  [x] Shear centre (Trefftz’s method)
   -  [x] Shear areas (global axis)
   -  [x] Shear areas (principal axis)
-  [x] Cross-section stresses

Solver:
~~~~~~~

-  [x] Direct solver
-  [x] CGS iterative solver
-  [x] Sparse matrices

Post-Processor:
~~~~~~~~~~~~~~~

- [x] Plot geometry
- [x] Plot mesh
- [x] Plot centroids
- [x] Plot cross-section stresses
- [ ] Generate cross-section report
- [ ] Export to Paraview

Additional Modules:
~~~~~~~~~~~~~~~~~~~

- [ ] Optimisation
- [ ] Reinforced Concrete
- [ ] Steel

.. _documentation: https://sectionproperties.readthedocs.io/
.. _installation wheel: https://www.lfd.uci.edu/~gohlke/pythonlibs/#meshpy

.. |Build Status| image:: https://travis-ci.com/robbievanleeuwen/section-properties.svg?branch=master
   :target: https://travis-ci.com/robbievanleeuwen/section-properties
.. |Documentation Status| image:: https://readthedocs.org/projects/sectionproperties/badge/?version=latest
   :target: https://sectionproperties.readthedocs.io/en/latest/?badge=latest
