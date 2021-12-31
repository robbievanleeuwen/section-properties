"""
Examples of using bridgebeams - using Spyder cells
"""

import sectionproperties.pre.library.bridge_sections as bridgebeams
import sectionproperties.pre.sections as sections
from sectionproperties.pre.pre import Material
from sectionproperties.analysis.cross_section import CrossSection

#%% Super-T

Dslab, w, t_f = 180, 2100, 75
super_t = bridgebeams.SuperTGirderSection(girder_type=5, w=w)
slab = sections.RectangularSection(Dslab, w, shift=[-w / 2, t_f])

precast = Material(
    name="65Mpa",
    elastic_modulus=37.4e3,
    poissons_ratio=0.2,
    yield_strength=65,
    density=2.4e-6,
    color="grey",
)
insitu = Material(
    name="40Mpa",
    elastic_modulus=32.8e3,
    poissons_ratio=0.2,
    yield_strength=40,
    density=2.4e-6,
    color="lightgrey",
)

geometry = sections.MergedSection([super_t, slab])
geometry.add_hole([0, -Dslab])
geometry.clean_geometry()
geometry.plot_geometry()

mesh = geometry.create_mesh(mesh_sizes=[500.0, 500.0])
section = CrossSection(geometry, mesh, materials=[precast, insitu])
section.plot_mesh(materials=True, alpha=0.4)

section.calculate_geometric_properties()
section.calculate_warping_properties()
section.display_results(fmt=".3f")

#%% I Girder - check table in AS5100.5 Fig. D1(A)

import pandas as pd

df = pd.DataFrame(columns=["Ag", "Zt", "Zb", "I", "dy", "th"])

for i in range(4):
    geometry = bridgebeams.IGirderSection(girder_type=i + 1)
    dims = geometry.get_girder_dims(girder_type=i + 1)
    d = sum(dims[-5:])
    mesh = geometry.create_mesh(mesh_sizes=[200.0])
    section = CrossSection(geometry, mesh)
    section.calculate_geometric_properties()
    section.calculate_warping_properties()

    A = section.get_area()
    th = A / (section.get_perimeter() / 2)

    df.loc[i] = [
        A,
        *(section.get_z()[:2]),
        section.get_ic()[0],
        d + section.get_c()[1],
        th,
    ]

print(df)


#%% I Girder and slab, modular ratio

Dslab, w = 175, 1400
i_girder = bridgebeams.IGirderSection(girder_type=4)
slab = sections.RectangularSection(Dslab, w, shift=[-w / 2, 0])

precast = Material(
    name="40Mpa",
    elastic_modulus=1,  # 32.8e3,
    poissons_ratio=0.2,
    yield_strength=40,
    density=2.4e-6,
    color="grey",
)
insitu = Material(
    name="25Mpa",
    elastic_modulus=0.7,  # 26.7e3,
    poissons_ratio=0.2,
    yield_strength=25,
    density=2.4e-6,
    color="lightgrey",
)

geometry = sections.MergedSection([i_girder, slab])
geometry.clean_geometry()
geometry.plot_geometry()

mesh = geometry.create_mesh(mesh_sizes=[500.0, 500.0])
section = CrossSection(geometry, mesh, materials=[precast, insitu])
section.plot_mesh(materials=True, alpha=0.4)

section.calculate_geometric_properties()
section.calculate_warping_properties()
section.display_results(fmt=".3f")

print("Compare with AS5100.5 Tb. D4(A) - ds = 150 mm; alpha_c = 0.7")
print(
    f"Torsion constant: {section.get_j() / precast.elastic_modulus / 1e6:.0f}x10^6 mm3"
)

#%% I Girder Only

geometry = bridgebeams.IGirderSection(girder_type=4)

mesh = geometry.create_mesh(mesh_sizes=[500.0])
section = CrossSection(geometry, mesh)
section.plot_mesh()

section.calculate_geometric_properties()
section.calculate_warping_properties()
section.display_results(fmt=".3f")

print("Compare with AS5100.5 Tb. D4(A) - girder only")
print(f"Torsion constant: {section.get_j() / 1e6:.0f}x10^6 mm3")
