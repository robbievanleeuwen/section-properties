.. _label-advanced_geom:

Advanced Geometry creation
==========================

The below tutorial was created to demonstrate the creation of valid geometries
for section analysis by combining multiple shapes.

Some key points to remember:

1. Geometries of two _different_ materials should not overlap (can create un-predictable results)
2. If two geometries of the _same_ materials are overlapping, then you should perform a union on the two sections
3. Two different section geometries that share a common edge (facet) should also share the same nodes (do not leave "floating" nodes along common edges)

These are general points to remember for any finite element analysis.

Note: `sectionproperties` will not prevent the creation of these ambiguous sections. The flexibility of the new
pre-processing engine (shapely) allows for a wide variety of intermediate modelling steps but the user must ensure
that the final model is one that is appropriate for analysis.

Creating Merged sections
------------------------

For this example, we will create a custom section out of two similar "I" sections::
    
    import sectionproperties.pre.sections as sections
    import sectionproperties.analysis.cross_section as cross_section

    i_sec1 = sections.i_section(d=250, b=150, t_f=13, t_w=10, r=12, n_r=12)
    i_sec2 = i_sec1.rotate_section(45)

..  figure:: ../images/examples/i_sec1.png
    :align: center
    :scale: 75 %

    I-section 1 SVG representation

..  figure:: ../images/examples/i_sec2.png
    :align: center
    :scale: 75 %

    I-section 2 SVG representation

Assign a unique material to each geometry::

    from sectionproperties.pre.pre import Material

    mat1 = Material("Material_1", 200e3, 0.3, 400, "red")
    mat2 = Material("Material_2", 150e3, 0.2, 200, "blue") # Just some differing properties

    i_sec1.material = mat1
    i_sec2.material = mat2

Now, we can use the `+` operator to naively combine these two sections into a `CompoundGeometry`::

    i_sec1 + i_sec2

..  figure:: ../images/examples/basic_compound.png
    :align: center
    :scale: 75 %

    SVG representation of our combined sections, note the two different materials.

However, if we use `.plot_geometry()` we will see that, even though we have two materials, we only have one control point for both geometries::

    (i_sec1 + i_sec2).plot_geometry()

..  figure:: ../images/examples/basic_compound_plot.png
    :align: center
    :scale: 75 %

    A naive compound of two geometries, each with different materials, but with only one analytical region.

If we went a few steps further and created a mesh and then plotted that mesh as part of an analysis section, we would see the unpredictable result of the mesh::

    cross_section.Section((i_sec1 + i_sec2).create_mesh([10])).plot_mesh()

..  figure:: ../images/examples/basic_combined_mesh_error.png
    :align: center
    :scale: 75 %

    The material information in the mesh did not turn out well!

To prevent ambiguity, there are a few options we can take. We can perform a simple union operation but that will lose
the material information for one of our sections, whichever section comes first in the operation::

..  figure:: ../images/examples/basic_union.png
    :align: center
    :scale: 75 %

    Combining two sections with union. Note how we only have `Material_2` now because `i_sec2` took precedence in the operation.

However, this is unsatisfactory as a solution. We want this section to more aptly represent a real section that might be created by cutting and welding two sections together.

Lets say we want the upright "I" section to be our main section and the diagonal section will be added on to it. 
The basic approach is require set operations performed in a few steps.

The quick way to do this would be use a difference operation and then combine it to the main section::

    (i_sec2 - i_sec1) + i_sec1

..  figure:: ../images/examples/combined_section_lucky.png
    :align: center
    :scale: 75 %

    Combining the sections this way _appears_ to give the result we want. However, this is a "lucky" combination.

Here is the plot of this section. You can now see five distinct regions demarcated by the five control points that were generated::

..  figure:: ../images/examples/combined_section_lucky_plot.png
    :align: center
    :scale: 75 %
    
However, this is a "lucky" combination. It's lucky because the regions where sections 1 and 2 share and edge,
they do not have nodes in common: the intersection nodes only exist on section 2 and not on section 1
(which still just contains its original nodes).

To combine these sections so that all intersection nodes are held in common, an extra step is required::

    cut_2_from_1 = (i_sec1 - i_sec2) # locates intersection nodes
    sec_1_nodes_added = cut_2_from_1 | sec_1

    # This can also be done in one line
    sec_1_nodes_added = (i_sec1 - i_sec2) | i_sec1

Now, when we use `.plot_geometry()`, we can see the additional nodes added to section 1::

    sec_1_nodes_added.plot_geometry()

..  figure:: ../images/examples/sec1_nodes_added.png
    :align: center
    :scale: 75 %

    The additional nodes from the cut portion are now merged as part of the section 1 geometry.

At this point, we can use our "section 1 with additional nodes" to create our complete geometry::

    analysis_geom = (i_sec2 - i_sec1) + sec_1_nodes_added
    analysis_geom.create_mesh([10])
    analysis_sec = cross_section.Section(analysis_geom)
    analysis_sec.plot_mesh()

..  figure:: ../images/examples/complete_combined_mesh.png
    :align: center
    :scale: 75 %

    The completed mesh and analysis section


