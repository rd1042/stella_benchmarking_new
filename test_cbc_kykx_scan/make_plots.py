""" """

import sys
sys.path.append("../postprocessing_tools")
from plotting_helper import make_comparison_plots, plot_gmvus, plot_gzvs, make_ky_scan_plots
from helper_ncdf import view_ncdf_variables, extract_data_from_ncdf
import matplotlib.pyplot as plt
import numpy as np
import glob
import re

# stella
stella_folder = "stella_kykx_scan"
stella_me1_folder = "stella_me1_kykx_scan"
stella_adiabatic_folder = "stella_adiabatic_kykx_scan"


# gs2
gs2_folder = "gs2_kykx_scan"
gs2_me1_folder = "gs2_me1_kykx_scan"
gs2_adiabatic_folder = "gs2_adiabatic_kykx_scan"



def get_sim_longnames(folder_longname, kx_str):
    """ """
    # Find all the input files in the folder.
    # NB might be more failsafe to find all sims with a .in, .out and
    # .omega file, but this is more work.
    sim_infile_longnames = glob.glob(folder_longname + "/*" + kx_str + ".in")
    unsorted_longnames = []
    ky_vals = []
    for sim_idx, sim_infile_longname in enumerate(sim_infile_longnames):
        # Get the sim longanme and ky.
        sim_longname = re.split(".in", sim_infile_longname)[0]
        unsorted_longnames.append(sim_longname)
        kymin_kymax_str = re.split("ky_kx_", sim_longname)[-1]
        ky_str = re.split("_", kymin_kymax_str)[0]
        ky_vals.append(float(ky_str))

    # Sort into ascending order of ky
    ky_vals = np.array(ky_vals)
    sort_idxs = np.argsort(ky_vals)
    ky_vals = ky_vals[sort_idxs]
    sim_longnames = []
    for sort_idx in sort_idxs:
        sim_longnames.append(unsorted_longnames[sort_idx])

    return ky_vals, sim_longnames

def compare_ky_scans(stella_folder, gs2_folder, kx_str_stella, kx_str_gs2, save_name, plot_apar=False, plot_bpar=False):
    """ """
    stella_ky, stella_longnames = get_sim_longnames(stella_folder, kx_str_stella)
    gs2_ky, gs2_longnames = get_sim_longnames(gs2_folder, kx_str_gs2)

    print("ky = ", stella_ky)
    print("stella_longnames = ", stella_longnames)
    print("gs2_longnames = ", gs2_longnames)

    # Expect gs2_ky = stella_ky; stop if not (to implement: do something
    # clever in this case)
    if np.max(abs(stella_ky - gs2_ky)) > 1e-4:
        print("Error! GS2 ky != stella ky . Stopping")
        print("stella_ky = ", stella_ky)
        print("gs2_ky = ", gs2_ky)
        sys.exit()
    make_ky_scan_plots(stella_longnames, gs2_longnames, stella_ky, save_name,
            gs2_pickle=None,  plot_apar=False, plot_bpar=False, plot_format=".png")

def plot_different_ky_scans():
    """ """

    # compare_ky_scans(stella_fapar1_fbpar0_me1_folder,
    #                    gs2_fapar1_fbpar0_me1_folder,
    #                    "images/fapar1_fbpar0")
    compare_ky_scans(stella_folder,
                       gs2_folder,
                       "_0.0_0.0", "_0.0",
                       "images/gyro_electrons_kx=0")
    # compare_ky_scans(stella_me1_folder,
    #                    gs2_me1_folder,
    #                    "_0.0_0.0", "_0.0",
    #                    "images/m1_electrons_kx=0")
    # compare_ky_scans(stella_adiabatic_folder,
    #                    gs2_adiabatic_folder,
    #                    "_0.0_0.0", "_0.0",
    #                    "images/adiabatic_electrons_kx=0")


    return

if __name__ == "__main__":
    plot_different_ky_scans()
