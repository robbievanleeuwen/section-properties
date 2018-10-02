import sectionproperties.pre.sections as sections
import sectionproperties.pre.pre as pre
from sectionproperties.analysis.cross_section import CrossSection

# create cross-section geometry
pfc_left = sections.PfcSection(d=200, b=75, t_f=12, t_w=6, r=12, n_r=8)
pfc_right = sections.PfcSection(d=200, b=75, t_f=12, t_w=6, r=12, n_r=8)
pfc_right.mirror_section(axis='y', mirror_point=[75, 0])
plate_top = sections.RectangularSection(d=16, b=250, shift=[-50, 200])
geometry = sections.MergedSection([pfc_left, pfc_right, plate_top])
geometry.add_hole([75, 100])
geometry = pre.GeometryCleaner(geometry).clean_geometry()
geometry.plot_geometry()

# create cross-section mesh
mesh = geometry.create_mesh(mesh_sizes=[5, 5, 10])

# create cross-section object
section = CrossSection(geometry, mesh)
section.display_mesh_info()
section.plot_mesh()

# calculate cross-section properties
section.calculate_geometric_properties(time_info=True)
section.calculate_warping_properties(time_info=True)
section.display_results()
section.plot_centroids()

# calculate cross-section stresses
stress_result = section.calculate_stress(Mxx=-20e6, Vy=-10e3, Mzz=5e6,
                                         time_info=True)

# show normal stresses
stress_result.plot_stress_zz()

# show torsion stresses
stress_result.plot_stress_mzz_zxy()
stress_result.plot_vector_mzz_zxy()

# show shear stresses
stress_result.plot_stress_v_zxy()
stress_result.plot_vector_v_zxy()

# show von-mises stress
stress_result.plot_stress_vm()
