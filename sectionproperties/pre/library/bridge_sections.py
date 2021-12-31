"""
This is a port of the bridgebeams package that was external to v1 sectionproperties
"""

import numpy as np
from sectionproperties.pre.geometry import Geometry


class SuperTGirderSection(Geometry):
    """Constructs a Super T Girder section to AS5100.5

    :param int girder_type: Type of Super T (1 to 5)
    :param int girder_subtype: Era Super T (1: pre-2001, 2:contemporary)
    :param float w: Overall width of top flange
    :param float t_w: Web thickness of the Super-T section (defaults to those of AS5100.5 Tb D3(B))
    :param float t_f: Thickness of top flange (VIC (default) = 75 mm; NSW = 90 mm)
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a T5 Super-T section with a 180 overlay slab and takes acount of
    the different material properties::

        import bridgebeams
        import sectionproperties.pre.sections as sections
        from sectionproperties.pre.pre import Material
        from sectionproperties.analysis.cross_section import CrossSection

        Dslab, w, t_f = 180, 2100, 75
        super_t = bridgebeams.SuperTGirderSection(girder_type=5, w=w)
        slab = sections.RectangularSection(Dslab, w, shift=[-w / 2, t_f])

        precast = Material(
            name="65Mpa",
            elastic_modulus=37.4e3,
            poissons_ratio=0.2,
            yield_strength=65,
            color="grey",
        )
        insitu = Material(
            name="40Mpa",
            elastic_modulus=32.8e3,
            poissons_ratio=0.2,
            yield_strength=40,
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

    Note that the properties are reported as ``modulus weighted'' properties (e.g. E.A) and can
    be normalized to the reference material by dividing by that elastic modulus::

        A_65 = section.get_ea() / precast.elastic_modulus

    The reported section centroids are already weighted.

    """

    def __init__(
        self, girder_type, girder_subtype=2, w=2100, t_w=None, t_f=75, shift=(0, 0)
    ):
        """Inits the SuperTGirderSection class"""

        if girder_type < 1 or girder_type > 5:
            msg = "Super-T Girder Type must be between 1 and 5"
            raise Exception(msg)

        if girder_subtype not in [1, 2]:
            msg = "Only 2 subtypes: pre- and post-2001"
            raise Exception(msg)

        if girder_type == 5 and girder_subtype == 1:
            msg = "Only girders T1 to T4 before 2001"
            raise Exception(msg)

        # Some super-t constants, refer to AS5100.5, Figure D1(B)
        d_fillet, b_fillet = 75, 100
        d_recess, b_recess = 25, 25
        d_chamfer, b_chamfer = 13, 13
        web_slope = 10.556
        flg_slope = 5.347
        w_nom = 1027
        if girder_subtype == 1:
            flg_slope = 5.0
            w_nom = 920

        # Dims for the specific girder type
        d, t_b, t_w_nom = self.get_girder_dims(girder_type)
        # Overriding default web thickness?
        if t_w is None:
            t_w = t_w_nom

        # assign control point as middle of bottom flange
        control_points = [[0, -d + t_b / 2]]
        super().__init__(control_points, shift)

        # Origin is middle at level of bottom of top flange
        # Some geometrics of the slope
        web_hyp = np.sqrt(1 + web_slope ** 2)
        web_horiz = t_w * web_slope / web_hyp
        x_fillet = w_nom / 2 - d_fillet * 1 / web_hyp

        # Right recess
        pt_b = w_nom / 2 - web_horiz + (t_f - d_recess) / web_slope
        pt_a = pt_b + b_recess
        self.points.append([pt_b, t_f - d_recess])
        self.points.append([pt_a, t_f - d_recess])
        self.points.append([pt_a, t_f])

        # Right flange
        self.points.append([w / 2, t_f])
        self.points.append([w / 2, d_chamfer])
        self.points.append([w / 2 - b_chamfer, 0])
        self.points.append([w_nom / 2 + b_fillet, 0])
        self.points.append([x_fillet, -d_fillet])

        # Bottom outer
        btm_corner = x_fillet - d / web_slope
        self.points.append([btm_corner + d_chamfer / web_slope, -(d - d_chamfer)])
        self.points.append([btm_corner - b_chamfer, -d])
        self.points.append([-btm_corner + b_chamfer, -d])
        self.points.append([-btm_corner - d_chamfer / web_slope, -(d - d_chamfer)])

        # Left flange
        self.points.append([-x_fillet, -d_fillet])
        self.points.append([-w_nom / 2 - b_fillet, 0])
        self.points.append([-w / 2 + b_chamfer, 0])
        self.points.append([-w / 2, d_chamfer])
        self.points.append([-w / 2, t_f])

        # Left Recess
        self.points.append([-pt_a, t_f])
        self.points.append([-pt_a, t_f - d_recess])
        self.points.append([-pt_b, t_f - d_recess])

        # Bottom inner - find intersection point
        y10 = t_f - d_recess
        y20 = -(d - t_b)
        y_inner = (1 / (1 / web_slope - flg_slope)) * (
            y10 / web_slope - y20 * flg_slope - pt_b
        )
        x_inner = flg_slope * (y20 - y_inner)
        self.points.append([x_inner, y_inner])
        self.points.append([0, y20])
        self.points.append([-x_inner, y_inner])

        # build facet list
        num_points = int(len(self.points))
        for i in range(num_points):
            # if we are not at the last point
            if i != num_points - 1:
                self.facets.append([i, i + 1])
            # if we are at the last point, complete the loop
            else:
                self.facets.append([i, 0])

        self.perimeter = list(range(len(self.facets)))

        self.shift_section()

    def get_girder_dims(self, girder_type):
        """Returns a dictionary of Super-T dimensions
        refer to AS5100.5, Appendix D
        """

        girder_dims = {
            "T1": {"d": 675, "t_b": 240, "t_w": 100},
            "T2": {"d": 925, "t_b": 240, "t_w": 100},
            "T3": {"d": 1125, "t_b": 260, "t_w": 100},
            "T4": {"d": 1425, "t_b": 260, "t_w": 100},
            "T5": {"d": 1725, "t_b": 325, "t_w": 120},
        }

        key = f"T{girder_type}"

        # Rather delicious code that assigns the values to the keys as variables
        for key, val in girder_dims.items():
            exec(key + "=val")

        d = girder_dims[key]["d"]
        t_b = girder_dims[key]["t_b"]
        t_w = girder_dims[key]["t_w"]

        return d, t_b, t_w


class IGirderSection(Geometry):
    """Constructs a precast I girder section to AS5100.5.

    As an example, replicate the table shown in AS5100.5 Fig. D1(A)::

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

    Note that the section depth is obtained by summing the heights from the
    section dictionary in `get_girder_dims`.

    """

    def __init__(self, girder_type, shift=(0, 0)):
        """Inits the SuperTGirderSection class"""

        if girder_type < 1 or girder_type > 4:
            msg = "I Girder Type must be between 1 and 4"
            raise Exception(msg)

        b_tf, b_bf, b_w, h_tf, h_ts, h_w, h_bs, h_bf = self.get_girder_dims(girder_type)

        # Some section constants
        d = sum([h_tf, h_ts, h_w, h_bs, h_bf])
        inset_tf = (b_tf - b_w) / 2
        inset_bf = (b_bf - b_w) / 2

        # assign control point as middle of bottom flange
        control_points = [[0, -d + h_bf / 2]]
        super().__init__(control_points, shift)

        # Origin at centre top; clockwise from top right corner
        self.points.append([b_tf / 2, 0])
        self.points.append([b_tf / 2, -h_tf])
        self.points.append([b_tf / 2 - inset_tf, -h_tf - h_ts])
        self.points.append([b_tf / 2 - inset_tf, -h_tf - h_ts - h_w])
        self.points.append([b_bf / 2, -d + h_bf])
        self.points.append([b_bf / 2, -d])
        self.points.append([-b_bf / 2, -d])
        self.points.append([-b_bf / 2, -d + h_bf])
        self.points.append([-b_tf / 2 + inset_tf, -h_tf - h_ts - h_w])
        self.points.append([-b_tf / 2 + inset_tf, -h_tf - h_ts])
        self.points.append([-b_tf / 2, -h_tf])
        self.points.append([-b_tf / 2, 0])

        # build facet list
        num_points = int(len(self.points))
        for i in range(num_points):
            if i != num_points - 1:  # not last point
                self.facets.append([i, i + 1])
            else:  # last point, so close
                self.facets.append([i, 0])

        self.perimeter = list(range(len(self.facets)))

        self.shift_section()

    def get_girder_dims(self, girder_type):
        """Returns a dictionary of I girder dimensions
        refer to AS5100.5, Appendix D
        """

        girder_dims = {
            "T1": {
                "b_tf": 200,
                "b_bf": 300,
                "b_w": 120,
                "h_tf": 100,
                "h_ts": 40,
                "h_w": 420,
                "h_bs": 90,
                "h_bf": 100,
            },
            "T2": {
                "b_tf": 350,
                "b_bf": 450,
                "b_w": 150,
                "h_tf": 100,
                "h_ts": 100,
                "h_w": 450,
                "h_bs": 150,
                "h_bf": 100,
            },
            "T3": {
                "b_tf": 450,
                "b_bf": 500,
                "b_w": 150,
                "h_tf": 130,
                "h_ts": 150,
                "h_w": 545,
                "h_bs": 175,
                "h_bf": 150,
            },
            "T4": {
                "b_tf": 500,
                "b_bf": 650,
                "b_w": 150,
                "h_tf": 150,
                "h_ts": 175,
                "h_w": 650,
                "h_bs": 250,
                "h_bf": 175,
            },
        }

        key = f"T{girder_type}"

        b_tf = girder_dims[key]["b_tf"]
        b_bf = girder_dims[key]["b_bf"]
        b_w = girder_dims[key]["b_w"]
        h_tf = girder_dims[key]["h_tf"]
        h_ts = girder_dims[key]["h_ts"]
        h_w = girder_dims[key]["h_w"]
        h_bs = girder_dims[key]["h_bs"]
        h_bf = girder_dims[key]["h_bf"]

        return b_tf, b_bf, b_w, h_tf, h_ts, h_w, h_bs, h_bf
