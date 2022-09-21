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

folder_name_dkh_explicit = "sims/wout_DKH_0_norm_beta_scan_explicit"

if __name__ == "__main__":
    print("Hello world")
    beta_vals = np.linspace(0,0.03, 10)
    # ky_vals = [0.05, 0.5, 1, 1.5, 3, 4.5, 6]
    make_beta_scan(folder_name_dkh_explicit, fprim_vals, tprim_vals, ky_vals)
