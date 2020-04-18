import numpy as np


def validate_properties(test, validation_list, section):
    """a"""

    for entry in validation_list:
        if entry["tol"] is None:
            tol = 1e-5
        else:
            tol = entry["tol"]

        if entry["prop"] == "area":
            val = section.get_area()
        elif entry["prop"] == "perimeter":
            val = section.get_perimeter()
        elif entry["prop"] == "ea":
            val = section.get_ea()
        elif entry["prop"] == "qx":
            (val, _) = section.get_q()
        elif entry["prop"] == "qy":
            (_, val) = section.get_q()
        elif entry["prop"] == "cx":
            (val, _) = section.get_c()
        elif entry["prop"] == "cy":
            (_, val) = section.get_c()
        elif entry["prop"] == "ixx_g":
            (val, _, _) = section.get_ig()
        elif entry["prop"] == "iyy_g":
            (_, val, _) = section.get_ig()
        elif entry["prop"] == "ixy_g":
            (_, _, val) = section.get_ig()
        elif entry["prop"] == "ixx_c":
            (val, _, _) = section.get_ic()
        elif entry["prop"] == "iyy_c":
            (_, val, _) = section.get_ic()
        elif entry["prop"] == "ixy_c":
            (_, _, val) = section.get_ic()
        elif entry["prop"] == "zxx_plus":
            (val, _, _, _) = section.get_z()
        elif entry["prop"] == "zxx_minus":
            (_, val, _, _) = section.get_z()
        elif entry["prop"] == "zyy_plus":
            (_, _, val, _) = section.get_z()
        elif entry["prop"] == "zyy_minus":
            (_, _, _, val) = section.get_z()
        elif entry["prop"] == "rx":
            (val, _) = section.get_rc()
        elif entry["prop"] == "ry":
            (_, val) = section.get_rc()
        elif entry["prop"] == "phi":
            val = section.get_phi()
        elif entry["prop"] == "i11_c":
            (val, _) = section.get_ip()
        elif entry["prop"] == "i22_c":
            (_, val) = section.get_ip()
        elif entry["prop"] == "z11_plus":
            (val, _, _, _) = section.get_zp()
        elif entry["prop"] == "z11_minus":
            (_, val, _, _) = section.get_zp()
        elif entry["prop"] == "z22_plus":
            (_, _, val, _) = section.get_zp()
        elif entry["prop"] == "z22_minus":
            (_, _, _, val) = section.get_zp()
        elif entry["prop"] == "r11":
            (val, _) = section.get_rp()
        elif entry["prop"] == "r22":
            (_, val) = section.get_rp()
        elif entry["prop"] == "x_pc":
            (val, _) = section.get_pc()
        elif entry["prop"] == "y_pc":
            (_, val) = section.get_pc()
        elif entry["prop"] == "x11_pc":
            (val, _) = section.get_pc_p()
        elif entry["prop"] == "y22_pc":
            (_, val) = section.get_pc_p()
        elif entry["prop"] == "sxx":
            (val, _) = section.get_s()
        elif entry["prop"] == "syy":
            (_, val) = section.get_s()
        elif entry["prop"] == "s11":
            (val, _) = section.get_sp()
        elif entry["prop"] == "s22":
            (_, val) = section.get_sp()
        elif entry["prop"] == "sf_xx_plus":
            (val, _, _, _) = section.get_sf()
        elif entry["prop"] == "sf_xx_minus":
            (_, val, _, _) = section.get_sf()
        elif entry["prop"] == "sf_yy_plus":
            (_, _, val, _) = section.get_sf()
        elif entry["prop"] == "sf_yy_minus":
            (_, _, _, val) = section.get_sf()
        elif entry["prop"] == "sf_11_plus":
            (val, _, _, _) = section.get_sf_p()
        elif entry["prop"] == "sf_11_minus":
            (_, val, _, _) = section.get_sf_p()
        elif entry["prop"] == "sf_22_plus":
            (_, _, val, _) = section.get_sf_p()
        elif entry["prop"] == "sf_22_minus":
            (_, _, _, val) = section.get_sf_p()
        elif entry["prop"] == "j":
            val = section.get_j()
        elif entry["prop"] == "gamma":
            val = section.get_gamma()
        elif entry["prop"] == "A_s11":
            (val, _) = section.get_As_p()
        elif entry["prop"] == "A_s22":
            (_, val) = section.get_As_p()
        elif entry["prop"] == "x11_se":
            (val, _) = section.get_sc_p()
        elif entry["prop"] == "y22_se":
            (_, val) = section.get_sc_p()
        else:
            raise KeyError("Incorrect property key: {0}".format(entry["prop"]))

        if entry["val"] != 0:
            calc_tol = (val - entry["val"]) / entry["val"]
        else:
            calc_tol = val - entry["val"]

        test.assertTrue(
            np.isclose(val, entry["val"], rtol=tol), msg="Prop: {0}; Tol: {1:.5e}".format(
                entry["prop"], calc_tol)
        )
