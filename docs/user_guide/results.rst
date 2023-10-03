Results
=======

Displaying the Results
----------------------

A list of cross-section properties that have been calculated by the performed analyses
can be printed to the terminal using the
:meth:`~sectionproperties.analysis.section.Section.display_results` method that belongs
to every  :class:`~sectionproperties.analysis.section.Section` object.

..  automethod:: sectionproperties.analysis.section.Section.display_results
    :noindex:

Retrieving Section Properties
-----------------------------

The best way to obtain the calculated cross-section properties is to use one of the
various ``get`` methods. As described in :ref:`label-material-affects-results`, the
results that can be retrieved depends on whether or not material properties have been
used in the analysis.

Non-Composite Analysis Results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following results can be obtained if no material properties have been applied, i.e.
only the *default material* has been used.

Geometric Analysis
""""""""""""""""""

.. autosummary::
    :nosignatures:

    ~sectionproperties.analysis.section.Section.get_area
    ~sectionproperties.analysis.section.Section.get_perimeter
    ~sectionproperties.analysis.section.Section.get_q
    ~sectionproperties.analysis.section.Section.get_ig
    ~sectionproperties.analysis.section.Section.get_c
    ~sectionproperties.analysis.section.Section.get_ic
    ~sectionproperties.analysis.section.Section.get_z
    ~sectionproperties.analysis.section.Section.get_rc
    ~sectionproperties.analysis.section.Section.get_ip
    ~sectionproperties.analysis.section.Section.get_phi
    ~sectionproperties.analysis.section.Section.get_zp
    ~sectionproperties.analysis.section.Section.get_rp

Warping Analysis
""""""""""""""""

.. autosummary::
    :nosignatures:

    ~sectionproperties.analysis.section.Section.get_j
    ~sectionproperties.analysis.section.Section.get_sc
    ~sectionproperties.analysis.section.Section.get_sc_p
    ~sectionproperties.analysis.section.Section.get_sc_t
    ~sectionproperties.analysis.section.Section.get_gamma
    ~sectionproperties.analysis.section.Section.get_as
    ~sectionproperties.analysis.section.Section.get_as_p
    ~sectionproperties.analysis.section.Section.get_beta
    ~sectionproperties.analysis.section.Section.get_beta_p

Plastic Analysis
""""""""""""""""

.. autosummary::
    :nosignatures:

    ~sectionproperties.analysis.section.Section.get_pc
    ~sectionproperties.analysis.section.Section.get_pc_p
    ~sectionproperties.analysis.section.Section.get_s
    ~sectionproperties.analysis.section.Section.get_sp
    ~sectionproperties.analysis.section.Section.get_sf
    ~sectionproperties.analysis.section.Section.get_sf_p

Composite Analysis Results
^^^^^^^^^^^^^^^^^^^^^^^^^^

The following results can be obtained if one or more material properties have been
applied. Some methods support obtaining transformed properties by supplying an optional
reference elastic modulus or material (typically those prefaced with an *e*, e.g.
``get_eic``). For an example of this see
`Retrieving Section Properties <../examples/results/get_results.ipynb>`_. Click on the
relevant ``get`` method for more information.

Geometric Analysis
""""""""""""""""""

.. autosummary::
    :nosignatures:

    ~sectionproperties.analysis.section.Section.get_area
    ~sectionproperties.analysis.section.Section.get_perimeter
    ~sectionproperties.analysis.section.Section.get_mass
    ~sectionproperties.analysis.section.Section.get_ea
    ~sectionproperties.analysis.section.Section.get_eq
    ~sectionproperties.analysis.section.Section.get_eig
    ~sectionproperties.analysis.section.Section.get_c
    ~sectionproperties.analysis.section.Section.get_eic
    ~sectionproperties.analysis.section.Section.get_ez
    ~sectionproperties.analysis.section.Section.get_rc
    ~sectionproperties.analysis.section.Section.get_eip
    ~sectionproperties.analysis.section.Section.get_phi
    ~sectionproperties.analysis.section.Section.get_ezp
    ~sectionproperties.analysis.section.Section.get_rp
    ~sectionproperties.analysis.section.Section.get_nu_eff
    ~sectionproperties.analysis.section.Section.get_e_eff
    ~sectionproperties.analysis.section.Section.get_g_eff

Warping Analysis
""""""""""""""""

.. autosummary::
    :nosignatures:

    ~sectionproperties.analysis.section.Section.get_ej
    ~sectionproperties.analysis.section.Section.get_sc
    ~sectionproperties.analysis.section.Section.get_sc_p
    ~sectionproperties.analysis.section.Section.get_sc_t
    ~sectionproperties.analysis.section.Section.get_egamma
    ~sectionproperties.analysis.section.Section.get_eas
    ~sectionproperties.analysis.section.Section.get_eas_p
    ~sectionproperties.analysis.section.Section.get_beta
    ~sectionproperties.analysis.section.Section.get_beta_p

Plastic Analysis
""""""""""""""""

.. autosummary::
    :nosignatures:

    ~sectionproperties.analysis.section.Section.get_pc
    ~sectionproperties.analysis.section.Section.get_pc_p
    ~sectionproperties.analysis.section.Section.get_mp
    ~sectionproperties.analysis.section.Section.get_mp_p

.. _label-material-affects-results:

How Material Properties Affect Results
--------------------------------------

``sectionproperties`` has been built in a generalised way to enable composite (multiple
material) analysis, see :ref:`label-theory-composite`. As a result, a number of
cross-section properties are calculated in a modulus-weighted manner. This means that
if materials are applied to geometries, it is assumed that the user is undertaking a
**composite** analysis and a number of the results can only be retrieved in a
modulus-weighted manner. If no materials are applied (i.e. the geometry has the
*default material*), it is assumed that the user is undertaking a **geometric-only**
(non-composite) analysis and the user can retrieve geometric section properties.

.. admonition:: Summary

    1. **Geometric-only analysis:** user does not provide any material properties,
       user can retrieve geometric section properties. For example, the user can use
       :meth:`~sectionproperties.analysis.section.Section.get_ic` to get the centroidal
       second moments of area.

    2. **Composite analysis:** user provides one or more material properties, user can
       retrieve geometric material property weighted properties. For example, the user
       can use :meth:`~sectionproperties.analysis.section.Section.get_eic` to get the
       modulus-weighted centroidal second moments of area.

    To illustrate this point, consider modelling a typical reinforced concrete section
    with a composite analysis approach. By providing material properties,
    ``sectionproperties`` calculates the gross section bending stiffness,
    :math:`(EI)_g`.

    .. math::
        (EI)_g = E_s \times I_s + E_c \times I_c

    This can be obtained using the
    :meth:`~sectionproperties.analysis.section.Section.get_eic` method. If the user
    wanted to obtain the transformed second moment of area for a code calculation, they
    could simply divide the gross bending stiffness by the elastic modulus for concrete:

    .. math::
        I_{c,eff} = \frac{(EI)_g}{E_c}

    This can be achieved in ``sectionproperties`` by providing an ``e_ref`` to
    :meth:`~sectionproperties.analysis.section.Section.get_eic`, for example:

    .. code-block:: python

        ei_gross = sec.get_eic()
        ic_eff = sec.get_eic(e_ref=concrete)

    For further detail, refer to the example in
    `Retrieving Section Properties <../examples/results/get_results.ipynb>`_.

.. admonition:: Modelling Recommendations
    :class: tip

    1. If there is only one material used in the geometry, *do not provide a material*
       and let ``sectionproperties`` use the default material.

    2. If there is only one material used in the geometry and the user is interested in
       material weighted properties, e.g. ``E.I`` or the plastic moment, *provide the
       material* and note that the results will be material property weighted.

    3. If there are multiple materials used in the geometry, *provide the materials* and
       note that the results will be material property weighted. If required, retrieve
       the cross-section properties using a reference elastic modulus or material to
       obtain transformed properties, which are often useful for design purposes.

Plotting Centroids
------------------

A plot of the various calculated centroids (i.e. elastic, plastic and shear centre), and
the principal axes can be produced by calling the
:meth:`~sectionproperties.analysis.section.Section.plot_centroids` method.

..  automethod:: sectionproperties.analysis.section.Section.plot_centroids
    :noindex:

Plotting Cross-Section Stresses
-------------------------------

After conducting a stress analysis of a cross-section based on applied actions, the
resulting stresses can be visualised using any of the
:meth:`~sectionproperties.post.stress_post.StressPost.plot_stress`,
:meth:`~sectionproperties.post.stress_post.StressPost.plot_stress_vector` or
:meth:`~sectionproperties.post.stress_post.StressPost.plot_mohrs_circles` methods.

Plot Stress Contours
^^^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.post.stress_post.StressPost.plot_stress
    :noindex:

.. admonition:: Example

    The following example plots a contour of the von Mises stress within a 150 x 90 x 12
    UA section resulting from the following actions:

    - :math:`N = 50` kN
    - :math:`M_{xx} = -5` kN.m
    - :math:`M_{22} = 2.5` kN.m
    - :math:`M_{zz} = 1.5` kN.m
    - :math:`V_x = 10` kN
    - :math:`V_y = 5` kN

    .. plot::
        :include-source: True
        :caption: Contour plot of the von Mises stress

        from sectionproperties.pre.library import angle_section
        from sectionproperties.analysis import Section

        # create geometry, mesh and section
        geom = angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
        geom.create_mesh(mesh_sizes=20)
        sec = Section(geometry=geom)

        # conduct analyses
        sec.calculate_geometric_properties()
        sec.calculate_warping_properties()
        stress = sec.calculate_stress(
            n=50e3, mxx=-5e6, m22=2.5e6, mzz=0.5e6, vx=10e3, vy=5e3
        )

        # plot stress contour
        stress.plot_stress(stress="vm", cmap="viridis", normalize=False)


Plot Stress Vectors
^^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.post.stress_post.StressPost.plot_stress_vector
    :noindex:

.. admonition:: Example

    The following example generates a vector plot of the shear stress within a 150 x 90
    x 12 UA section resulting from a torsion moment of 1 kN.m:

    .. plot::
        :include-source: True
        :caption: Contour plot of the von Mises stress

        from sectionproperties.pre.library import angle_section
        from sectionproperties.analysis import Section

        # create geometry, mesh and section
        geom = angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
        geom.create_mesh(mesh_sizes=20)
        sec = Section(geometry=geom)

        # conduct analyses
        sec.calculate_geometric_properties()
        sec.calculate_warping_properties()
        stress = sec.calculate_stress(mzz=1e6)

        # plot stress contour
        stress.plot_stress_vector(stress="mzz_zxy", cmap="viridis", normalize=False)

Plot Mohr's Circles
^^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.post.stress_post.StressPost.plot_mohrs_circles
    :noindex:

Retrieving Cross-Section Stresses
---------------------------------

Numerical values for cross-section stress can also be obtained with the
:meth:`~sectionproperties.analysis.section.Section.get_stress_at_points` and
:meth:`~sectionproperties.post.stress_post.StressPost.get_stress` methods.

Get Stress at Points
^^^^^^^^^^^^^^^^^^^^

This method can be used to obtain the stress at one or multiple points. A geometric
analysis must be performed prior to calling this method. Further, if the shear force or
torsion is non-zero, a warping analysis must also be performed. See
`Retrieving Stresses <../examples/results/get_stress.ipynb>`_ for an example of how this
method can be used to plot the stress distribution along a line.

..  automethod:: sectionproperties.analysis.section.Section.get_stress_at_points
    :noindex:

General Get Stress
^^^^^^^^^^^^^^^^^^

This method must be used after a stress analysis. It returns a data structure containing
all the stresses within the cross-section, grouped by material.

..  automethod:: sectionproperties.post.stress_post.StressPost.get_stress
    :noindex:
