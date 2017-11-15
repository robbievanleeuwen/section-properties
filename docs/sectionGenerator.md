# Cross-section generators

The following functions can be used to generate a standard cross-section. N.B. for fillet radii, if a radius of zero is desired, enter a radius of zero and number of points (n_r) equal to one. If you have a request for a helper function to be added, don't hesitate to contact me. Refer to the examples on my [blog](https://robbievanleeuwen.github.io) for more information.

### Combined Section Generator

```python
(points, facets, holes) = sectionGenerator.combineShapes(sections)
'''
Combines multiple sections (as long as there is connectivity between all
elements). Takes a list of dictionaries as an input. Each dictionary defines
a section and should have the following keys:
    - section1['points']: list of section nodes
    - section1['facets']: list of section facets
    - section1['holes']: list of section holes
    - section1['x']: x-offset for section
    - section1['y']: y-offset for section
'''
```

### Circular Hollow Section (CHS)

```python
(points, facets, holes) = sectionGenerator.CHS(d, t, n)
'''
Constructs a circular hollow section with diameter d, thickness t and uses
n points to construct the inner and outer circles.
'''
```

### Rectangular Hollow Section (RHS)

```python
(points, facets, holes) = sectionGenerator.RHS(d, b, t, r_out, n_r)
'''
Constructs a rectangular hollow section with depth d, width b, thickness t,
outer radius r_out, using n_r points to construct the inner and outer radii.
'''
```

### Rectangular Hollow Section with Split (RHS)

```python
(points, facets, holes) = sectionGenerator.RHS_Split(d, b, b_split, t, r_out, n_r)
'''
  Constructs a rectangular hollow section with depth d, width b, split
  thickness b_split, thickness t, outer radius r_out, using n_r points to
  construct the inner and outer radii.
  '''
```

### I-Section (e.g. UB or UC)

```python
(points, facets, holes) = sectionGenerator.ISection(d, b, tf, tw, r, n_r)
'''
Constructs an I-section with depth d, width b, flange thickness tf, web
thickness tw, root radius r, using n_r points to construct the root radius.
'''
```

### Channel Section (PFC)

```python
(points, facets, holes) = sectionGenerator.PFC(d, b, tf, tw, r, n_r)
'''
Constructs a PFC section with depth d, width b, flange thickness tf, web
thickness tw, root radius r, using n_r points to construct the root radius.
'''
```

### Angle Section (e.g. EA or UA)

```python
(points, facets, holes) = sectionGenerator.Angle(d, b, t, r_root, r_toe, n_r)
'''
Constructs an angle section with depth d, width b, thickness t, root radius
r_root, toe radius r_toe using n_r points to construct the root radius.
'''
```

### Rectangular Section

```python
(points, facets, holes) = sectionGenerator.Flat(d, b)
'''
Constructs a rectangular section with depth d and width b.
'''
```

### Circular Section

```python
(points, facets, holes) = sectionGenerator.Round(d, n)
'''
Constructs a solid cicular bar with diameter d, using n points to construct
the circle.
'''
```

### Tee Section (e.g. BT or CT)

```python
(points, facets, holes) = sectionGenerator.Tee(d, b, tf, tw, r, n_r)
'''
Constructs a Tee section with depth d, width b, flange thickness tf, web
thickness tw, root radius r, using n_r points to construct the root radius.
'''
```

### Cee Section

```python
(points, facets, holes) = sectionGenerator.Cee(d, b, l, t, r_out, n_r)
'''
Constructs a Cee section with depth d, width b, lip l, thickness t, outer
radius r_out, using n_r points to construct the root radius.
'''
```

### Zed Section

```python
(points, facets, holes) = sectionGenerator.Zed(d, b1, b2, l, t, r_out, n_r)
'''
Constructs a Zed section with depth d, left flange width b1, right flange
width b2, lip l, thickness t, outer radius r_out, using n_r points to
construct the root radius.
'''
```

### Cruciform Section

```python
(points, facets, holes) = sectionGenerator.Cruciform(d, b, t, r, n_r)
'''
Constructs a cruciform section with depth d, width b, thickness t, root
radius r, using n_r points to construct the root radius.
'''
```
