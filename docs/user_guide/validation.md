(label-testing)=

# Validation

``sectionproperties`` has a (slowly) growing suite of tests. The testing suite serves to
verify functionality and find exceptions from edge cases, but also validate the results
against known cases. Since this library performs engineering calculations, it should
have some manner of proving its accuracy. Each analyst who uses it is responsible for
their own projects and calculations, but having a high level of confidence that the
software can produce *correct* results, *given the correct inputs*, is a boon to all
users. Some test results and explanations from the latter category will be outlined on
this page, since the former really serves no use to the end user.

## Textbook Examples

An obvious starting location is replicating examples from academic texts.
"Aircraft Structures" by David J. Peery is a highly regarded text in the
field of aerospace structures [[1]](peery).

### Peery - Symmetric Sections

The simplest example of a realistic section problem is the symmetric I-Beam, with a
free-fixed boundary condition and a transverse tip load. The free-body-diagram and shear
and moment diagrams are shown in the problem statement, referenced below. This problem
is Example 1 in Section 6.2.

```{figure} ../_static/validation/peery_6-2-1.png
:width: 40%
:alt: From Peery, Sec. 6.2 Example 1
:align: center
```

#### Check #1

In ``sectionproperties``, there are multiple ways to set this problem up. We could
different shapely geometries and merge together, or a set of custom points, or a
built-in constructor. For the sake of simplicity, this simpler I-section is identical to
the Nastran I-section definition, so it makes sense to utilize the built-in constructor
from {func}`pre.library.nastran_sections.nastran_i()` (#TODO).

Using an arbitrarily coarse mesh, the properties can then be directly calculated from
the class method {meth}`Section.calculate_geomtric_properties()` (#TODO).

Peery lists the second moment of area about the primary bending axis as a value of 43.3
[in{math}`^4`]. For the automated tests in this library, we check against this hardcoded
value, with a tolerance of {math}`\pm 0.1`%.

#### Check #2

As a final check against this example, we can calculate the maximum bending stress on
the I-beam. From simple statics, the maximum moment from the FBD will be 800,000
[in-lbs] at the fixed end. Applying this moment to our section from before will allow
computation of stress over the FEM.

Peery quotes the peak value at 55.5 [ksi], which is rounded to the nearest decimal
place. From the equatio listed in the text, the theoretical value is actually 55,427.3
[psi].

```{math}
  f = \frac{My}{I} = 55,427.3 = 55,400
```

Again, the automated test against this checks the hardcoded value with a tolerance of
{math}`\pm 0.1`%. For accuracy, 55,427.3 is used instead of the rounded value.

For full details and the most updated code of this exampe, see the [examples page](#)
(#TODO) in the documentation gallery. For the exact test code execution, check the
[source](https://github.com/robbievanleeuwen/section-properties/blob/master/sectionproperties/tests/test_peery.py).

### Peery - Unsymmetric Sections

For a more complex example, we can turn to Example 1 in Section 7.2 of Peery. Here, we
have a still-simplified Z-section, but bending about two axes. Note axes definitions in
the problem statement. Beam axial direction in ``sectionproperties`` is always
referenced as the z-axis, and loads must be applied in this coordinate system.

```{figure} ../_static/validation/peery_7-2-1_1.png
:width: 40%
:alt: From Peery, Sec. 7.2 Example 1
:align: center
```

The construction of this geometry takes a similar approach to Ex 6.2.1, and utilizes a
built-in factory, {func}`pre.library.nastran_sections.nastran_zed()` (#TODO). The only
difference you may notice in the test code is usage of a custom class for ease of
initialization. This is not necessary.

Using an arbitrarily coarse mesh, the properties can then be directly calculated from
the class method {meth}`Section.calculate_geomtric_properties()` (#TODO). Each property
listed directly by Peery is taken as a hardcoded value and checked against, within the
testing suite.

```{eval-rst}
+---------------+--------------------------+
|   Property    |   Peery Value            |
+===============+==========================+
|  I_x          |   693.3 [in\ :sup:`4`]   |
+---------------+--------------------------+
|  I_y          |   173.3 [in\ :sup:`4`]   |
+---------------+--------------------------+
|  I_xy         |   -240  [in\ :sup:`4`]   |
+---------------+--------------------------+
|  I_p          |   787.1 [in\ :sup:`4`]   |
+---------------+--------------------------+
|  I_q          |    79.5 [in\ :sup:`4`]   |
+---------------+--------------------------+
|  theta        |   21.35 [deg]            |
+---------------+--------------------------+
```

For stress results, the theoretical values follow the biaxial bending equation. These
values are checked against automatically in the testing suite. Note that again Peery
rounds the values quoted directly, for simplicity. The testing suite also verifies that
the theoretical value as per the equation matches the theoretical value quoted in the
text, which also matches the computed value from the ``sectionproperties`` FEM.

```{math}}}
    f_b = \frac{M_x I_{xy} - M_y I_x}{I_x I_y - I_{xy}^2}x +
          \frac{M_y I_{xy} - M_x I_y}{I_x I_y - I_{xy}^2}y
```

```{eval-rst}
+---------+-----+-----+---------+---------+----------------------+
|  Point  |  x  |  y  |  -494x  |  -315y  |  :math:`f_b`, [psi]  |
+=========+=====+=====+=========+=========+======================+
|   A     | -5  |  4  |  2470   |  -1260  |  1210.0 = 1210       |
+---------+-----+-----+---------+---------+----------------------+
|   B     | -5  |  6  |  2470   |  -1890  |  580.0 = 580         |
+---------+-----+-----+---------+---------+----------------------+
|   C     |  1  |  6  |  -494   |  -1890  |  -2384.0 = -2380     |
+---------+-----+-----+---------+---------+----------------------+
```

```{figure} ../_static/validation/peery_7-2-1_2.png
:width: 40%
:alt: From Peery, Sec. 7.2 Example 1 results
:align: center
```

For full details and the most updated code of this exampe, see the [examples page](#)
(#TODO) in the documentation gallery. For the exact test code execution, check the
[source](https://github.com/robbievanleeuwen/section-properties/blob/master/sectionproperties/tests/test_peery.py).

## References

```{eval-rst}
.. _peery:
.. [1] D. J. Peery, *Aircraft Structures.* New York: Dover Publications, 2011.
    ISBN-10: 0-486-48580-3
```
