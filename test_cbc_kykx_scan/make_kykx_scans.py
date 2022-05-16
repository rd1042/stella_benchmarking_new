"""Code to generate scans in beta' for both stella (running on MARCONI)
and GS2 (running on VIKING) """

import numpy as np
import re
import sys
import os
import shutil

import sys
sys.path.append("../generate_sims")
from make_param_scans import construct_ky_kx_scan


## Define folder names here
# stella
stella_folder = "stella_kykx_scan"
stella_me1_folder = "stella_me1_kykx_scan"
stella_adiabatic_folder = "stella_adiabatic_kykx_scan"


# gs2
gs2_folder = "gs2_kykx_scan"
gs2_me1_folder = "gs2_me1_kykx_scan"
gs2_adiabatic_folder = "gs2_adiabatic_kykx_scan"


def make_kykx_scans():
    """Construct simulations scanning ky and kx.
    The test cases are:
    (1) Gyrokinetic electrons
    (2) Gyrokinetic electrons but with me=1
    (3) Adiabatic electrons """

    ky_vals1 = np.linspace(0, 5.0, 11)
    ky_vals2 = np.linspace(6.0, 10.0, 5)
    ky_vals3 = np.linspace(12, 40, 15)
    ky_vals = np.concatenate([ky_vals1, ky_vals2, ky_vals3])
    kx_vals = np.linspace(0.0, 1.0, 2)

    ##
    construct_ky_kx_scan(stella_folder, ky_vals, kx_vals, "stella")
    construct_ky_kx_scan(stella_me1_folder, ky_vals, kx_vals, "stella")
    construct_ky_kx_scan(stella_adiabatic_folder, ky_vals, kx_vals, "stella")
    # construct_ky_kx_scan(gs2_folder, ky_vals, kx_vals, "gs2")
    # construct_ky_kx_scan(gs2_me1_folder, ky_vals, kx_vals, "gs2")
    # construct_ky_kx_scan(gs2_adiabatic_folder, ky_vals, kx_vals, "gs2")

    return



if __name__ == "__main__":
    print("Hello world")
    make_kykx_scans()
