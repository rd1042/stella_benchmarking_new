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
from construct_stella_files_in_folders import make_beta_scan

folder_name_dkh_explicit_betalr = "sims/wout_DKH_0_norm_betalr_scan_explicit"
folder_name_dkh_implicit_betalr = "sims/wout_DKH_0_norm_betalr_scan_implicit"
folder_name_dkh_explicit_betalr_fbpar0 = "sims/wout_DKH_0_norm_betalr_scan_fbpar0_explicit"
folder_name_w003_explicit = "sims/w003_beta_scan_explicit"
folder_name_w003_implicit_ky1 = "sims/w003_beta_scan_implicit_ky1"
folder_name_w003_implicit_ky1_upw0 = "sims/w003_beta_scan_implicit_ky1_upw0"
folder_name_w003_implicit_ky2 = "sims/w003_beta_scan_implicit_ky2"
folder_name_w003_implicit_ky15 = "sims/w003_beta_scan_implicit_ky1.5"
folder_name_w003_explicit_fbpar0 = "sims/w003_beta_scan_explicit_fbpar0"
folder_name_w003_implicit_ky1_fbpar0 = "sims/w003_beta_scan_implicit_ky1_fbpar0"
folder_name_w003_implicit_ky2_fbpar0 = "sims/w003_beta_scan_implicit_ky2_fbpar0"
folder_name_w003_implicit_ky15_fbpar0 = "sims/w003_beta_scan_implicit_ky1.5_fbpar0"

if __name__ == "__main__":
    print("Hello world")
    beta_vals = np.linspace(0,0.1, 11)
    # make_beta_scan(folder_name_dkh_explicit_betalr, beta_vals)
    # make_beta_scan(folder_name_dkh_explicit_betalr_fbpar0, beta_vals)
    # make_beta_scan(folder_name_dkh_implicit_betalr, beta_vals)

    # make_beta_scan(folder_name_w003_explicit, beta_vals)
    # make_beta_scan(folder_name_w003_explicit_fbpar0, beta_vals)
    # make_beta_scan(folder_name_w003_implicit_ky1, beta_vals)
    # make_beta_scan(folder_name_w003_implicit_ky1_fbpar0, beta_vals)
    # make_beta_scan(folder_name_w003_implicit_ky15, beta_vals)
    # make_beta_scan(folder_name_w003_implicit_ky15_fbpar0, beta_vals)
    # make_beta_scan(folder_name_w003_implicit_ky2, beta_vals)
    # make_beta_scan(folder_name_w003_implicit_ky2_fbpar0, beta_vals)
    make_beta_scan(folder_name_w003_implicit_ky1_upw0, beta_vals)
