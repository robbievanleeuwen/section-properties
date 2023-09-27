.. _label-materials:

Materials
=========

Assigning materials to :class:`~sectionproperties.pre.geometry.Geometry` objects is
completely optional in ``sectionproperties``. In fact, if you are not conducting a
composite analysis it is recommended to not specify material properties as this adds
little value to the analysis results.

If undertaking a composite analysis, materials can be created using the
:class:`~sectionproperties.pre.pre.Material` object:

.. autoclass:: sectionproperties.pre.pre.Material
    :noindex:

:class:`~sectionproperties.pre.pre.Material` objects are assigned to
:class:`~sectionproperties.pre.geometry.Geometry` objects, learn more about how to
manipulate a :class:`~sectionproperties.pre.geometry.Geometry`'s material here,
:ref:`label-geom-material`.

Assigning materials affects the results reported by ``sectionproperties``, learn more
:ref:`here<label-material-affects-results>`.
