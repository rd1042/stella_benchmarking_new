"""Construct 3 sets of sims scanning fprim, tprim, ky:
fprim=linspace(0,9,8)
tprim=linspace(0,9,8)
ky = [0.05, 0.5, 1, 1.5, 3, 4.5, 6]
Do this for each equilibrium (003, 015, 027)
In each case, set ions and electrons equal."""

import numpy as np
import sys
import re
import math
from construct_stella_files_in_folders import make_fprim_tprim_ky_scan

folder_name_1 = "sims/w003_em_linear_fprim_tprim_explicit_scan"
folder_name_2 = "sims/w003_em_linear_fprim_tprim_implicit_scan"
folder_name_expl_higher_ky_em = "sims/w003_em_linear_fprim_tprim_explicit_scan_higher_ky"
folder_name_expl_higher_ky_es = "sims/w003_es_linear_fprim_tprim_explicit_scan_higher_ky"
folder_name_impl_higher_ky_em = "sims/w003_em_linear_fprim_tprim_implicit_scan_higher_ky"
folder_name_impl_higher_ky_es = "sims/w003_es_linear_fprim_tprim_implicit_scan_higher_ky"

if __name__ == "__main__":
    print("Hello world")
    fprim_vals = np.linspace(0,9,8)
    tprim_vals = np.linspace(0,9,8)
    # ky_vals = [0.05, 0.5, 1, 1.5, 3, 4.5, 6]
    # make_fprim_tprim_ky_scan(folder_name_1, fprim_vals, tprim_vals, ky_vals)
    # ky_vals = [0.5]
    # make_fprim_tprim_ky_scan(folder_name_1, fprim_vals, tprim_vals, ky_vals)
    # make_fprim_tprim_ky_scan(folder_name_2, fprim_vals, tprim_vals, ky_vals)
    ky_vals = [1.5, 3.0]
    # make_fprim_tprim_ky_scan(folder_name_expl_higher_ky_em, fprim_vals, tprim_vals, ky_vals)
    make_fprim_tprim_ky_scan(folder_name_expl_higher_ky_es, fprim_vals, tprim_vals, ky_vals)
    make_fprim_tprim_ky_scan(folder_name_impl_higher_ky_em, fprim_vals, tprim_vals, ky_vals)
    make_fprim_tprim_ky_scan(folder_name_impl_higher_ky_es, fprim_vals, tprim_vals, ky_vals)
