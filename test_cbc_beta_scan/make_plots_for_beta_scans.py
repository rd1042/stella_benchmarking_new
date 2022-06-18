""" """

import sys
sys.path.append("../postprocessing_tools")
from plotting_helper import make_comparison_plots, plot_gmvus, plot_gzvs, make_beta_scan_plots
from plotting_helper import make_beta_scan_plots_for_poster
from helper_ncdf import view_ncdf_variables, extract_data_from_ncdf
import matplotlib.pyplot as plt
import numpy as np
import glob
import re
import pickle

## Define folder names here
# stella
# stella_fapar0_fbpar1_me1_folder = "stella_fapar0_fbpar1_me1_beta_scan"
# stella_fapar1_fbpar0_me1_folder = "stella_fapar1_fbpar0_me1_beta_scan"
# stella_fapar1_fbpar1_me1_folder = "stella_fapar1_fbpar1_me1_beta_scan"
stella_fapar0_fbpar1_folder = "stella_fapar0_fbpar1_beta_scan"
stella_fapar1_fbpar0_folder = "stella_fapar1_fbpar0_beta_scan"
stella_fapar1_fbpar1_folder = "stella_fapar1_fbpar1_beta_scan"

# gs2
gs2_fapar0_fbpar1_me1_folder = "gs2_fapar0_fbpar1_me1_beta_scan"
gs2_fapar1_fbpar0_me1_folder = "gs2_fapar1_fbpar0_me1_beta_scan"
gs2_fapar1_fbpar1_me1_folder = "gs2_fapar1_fbpar1_me1_beta_scan"

stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt48_nvpa24_nmu18_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt48_nvpa24_nmu18_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt64_nvpa24_nmu18_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt64_nvpa24_nmu18_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt48_nvpa36_nmu18_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt48_nvpa36_nmu18_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt128_nvpa24_nmu18_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt128_nvpa24_nmu18_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_dt5em4_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_dt5em4_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_noupwind_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_noupwind_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_center_dgdz_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_center_dgdz_fapar1_fbpar1"


stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_ky_1_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_ky_1_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_ky_1_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_ky_1_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_ky_1_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_ky_1_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_ky_1_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_ky_1_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_ky_5_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_ky_5_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_ky_5_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_ky_5_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_ky_5_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_ky_5_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_ky_5_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_ky_5_fapar1_fbpar1"


gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_fapar1_fbpar1_folder="gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_fapar1_fbpar1"
gs2_beta_scan_ky_05_np3_nt48_ng8_ne18_fapar1_fbpar1_folder="gs2_beta_scan_ky_05_np3_nt48_ng8_ne18_fapar1_fbpar1"

gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_ky_1_fapar1_fbpar1_folder="gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_ky_1_fapar1_fbpar1"
gs2_beta_scan_ky_05_np3_nt48_ng8_ne18_ky_1_fapar1_fbpar1_folder="gs2_beta_scan_ky_05_np3_nt48_ng8_ne18_ky_1_fapar1_fbpar1"
gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_ky_5_fapar1_fbpar1_folder="gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_ky_5_fapar1_fbpar1"
gs2_beta_scan_ky_05_np3_nt48_ng8_ne18_ky_5_fapar1_fbpar1_folder="gs2_beta_scan_ky_05_np3_nt48_ng8_ne18_ky_5_fapar1_fbpar1"

compare_em_old = "old_em_stella_comparison/electromagnetic/"
compare_em_new = "old_em_stella_comparison/electromagnetic-new/"
compare_em_new_str_impl = "old_em_stella_comparison/electromagnetic-new_stream_implicit/"

old_stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_2_folder=compare_em_old + "stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_2"
old_stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2_folder=compare_em_old + "stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2"
old_stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2_folder=compare_em_old + "stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2"
old_stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2_folder=compare_em_old + "stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2"
new_stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_2_folder=compare_em_new + "stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_2"
new_stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2_folder=compare_em_new + "stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2"
new_stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2_folder=compare_em_new + "stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2"
new_stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2_folder=compare_em_new + "stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2"

stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2_streaming_implicit_folder = compare_em_new_str_impl + "stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2"
stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2_streaming_implicit_folder = compare_em_new_str_impl + "stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2"
stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2_streaming_implicit_folder = compare_em_new_str_impl + "stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2"
stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_streaming_implicit_folder = "stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_streaming_implicit"
stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_streaming_implicit_folder = "stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_streaming_implicit"
stella_beta_scan_ky_05_np2_nt64_nvpa24_nmu18_fapar1_fbpar1_streaming_implicit_folder = "stella_beta_scan_ky_05_np2_nt64_nvpa24_nmu18_fapar1_fbpar1_streaming_implicit"
stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1_streaming_implicit_folder = "stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1_streaming_implicit"
stella_beta_scan_bespoke_explicit_folder = "stella_bespoke_explicit"
stella_beta_scan_bespoke_stream_implicit_folder = "stella_bespoke_str_impl"

def get_sim_longnames(folder_longname):
    """ """
    # Find all the input files in the folder.
    # NB might be more failsafe to find all sims with a .in, .out and
    # .omega file, but this is more work.
    sim_infile_longnames = glob.glob(folder_longname + "/*.in")
    unsorted_longnames = []
    beta_vals = []
    for sim_idx, sim_infile_longname in enumerate(sim_infile_longnames):
        # Get the sim longanme and beta.
        sim_longname = re.split(".in", sim_infile_longname)[0]
        unsorted_longnames.append(sim_longname)
        beta_str = re.split("beta_", sim_longname)[-1]
        beta_vals.append(float(beta_str))

    # Sort into ascending order of beta
    beta_vals = np.array(beta_vals)
    sort_idxs = np.argsort(beta_vals)
    beta_vals = beta_vals[sort_idxs]
    sim_longnames = []
    for sort_idx in sort_idxs:
        sim_longnames.append(unsorted_longnames[sort_idx])

    return beta_vals, sim_longnames

def compare_beta_scans(stella_folder, gs2_folder, save_name, plot_apar=False, plot_bpar=False):
    """ """
    stella_beta, stella_longnames = get_sim_longnames(stella_folder)
    gs2_beta, gs2_longnames = get_sim_longnames(gs2_folder)

    print("stella_longnames = ", stella_longnames)
    print("gs2_longnames = ", gs2_longnames)

    # Expect gs2_beta = stella_beta; stop if not (to implement: do something
    # clever in this case)
    if np.max(abs(stella_beta - gs2_beta)) > 1e-4:
        print("Error! GS2 beta != stella beta . Stopping")
        print("stella_beta = ", stella_beta)
        print("gs2_beta = ", gs2_beta)
        sys.exit()
    make_beta_scan_plots(stella_longnames, gs2_longnames, gs2_beta, save_name,
            gs2_pickle=None,  plot_apar=True, plot_bpar=False, plot_format=".png")

def plot_different_beta_scans():
    """ """

    compare_beta_scans(stella_fapar1_fbpar0_me1_folder,
                       gs2_fapar1_fbpar0_me1_folder,
                       "images/fapar1_fbpar0")
    compare_beta_scans(stella_fapar1_fbpar1_me1_folder,
                       gs2_fapar1_fbpar1_me1_folder,
                       "images/fapar1_fbpar1")
    compare_beta_scans(stella_fapar0_fbpar1_me1_folder,
                       gs2_fapar0_fbpar1_me1_folder,
                       "images/fapar0_fbpar1")

    return

def plot_stella_scan_vs_gs2_pickle():
    """ """

    stella_beta, stella_longnames = get_sim_longnames(stella_fapar1_fbpar1_folder)
    # print("stella_longnames = ", stella_longnames)
    # print("stella_beta = ", stella_beta)
    # sys.exit()
    make_beta_scan_plots(stella_longnames, [], stella_beta, "images/beta_scan_fbpar1_fbpar1",
                        gs2_pickle="gs2_beta_scan/omega_values.pickle")
    return

def plot_stella_scan_vs_gs2_pickle_for_poster():
    """ """

    stella_beta, stella_longnames = get_sim_longnames(stella_fapar1_fbpar1_folder)
    # print("stella_longnames = ", stella_longnames)
    # print("stella_beta = ", stella_beta)
    # sys.exit()
    make_beta_scan_plots_for_poster(stella_longnames, stella_beta, "images/beta_scan_fbpar1_fbpar1_poster.png",
                        gs2_pickle="gs2_beta_scan/omega_values.pickle")
    return

def make_beta_plots_from_pickles(pickle_longnames, labels, marker_size=1, omega_diff=False):
    """Either plot Omega(beta), or, if omega_diff=True, plot (Omega-Omega_ref) (beta),
    where Omega_ref is taken from the first pickle. """

    marker_list = ["s", "o", "P", "X", "v", "^", "<", ">", "1", "2", "3"]
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)

    for folder_idx, pickle_longname in enumerate(pickle_longnames):
        myfile = open(pickle_longname, "rb")
        [beta_vals, gamma_vals, omega_vals] = pickle.load(myfile)
        myfile.close()
        if not omega_diff:
            ax1.plot(beta_vals, omega_vals, label=labels[folder_idx], marker=marker_list[folder_idx], mfc="none", markersize=marker_size)
            ax2.plot(beta_vals, gamma_vals, label=labels[folder_idx], marker=marker_list[folder_idx], mfc="none", markersize=marker_size)
        else:
            if folder_idx==0:
                # This is the reference beta
                beta_ref = beta_vals
                gamma_ref = gamma_vals
                omega_ref = omega_vals
            ## Whether or not this is the reference beta scan, plot Omega-Omega_ref
            # Find where beta vals match.
            beta_to_compare = []
            gamma_diff = []
            omega_diff = []
            for beta_idx, beta_val in enumerate(beta_vals):
                ## Find if the beta val is within a small tolerance of any
                ## beta_ref vals
                closest_beta_ref_idx = np.argmin(abs(beta_ref - beta_val))
                if abs(beta_ref[closest_beta_ref_idx] - beta_val) < 1e-6:
                    beta_to_compare.append(beta_val)
                    gamma_diff.append(gamma_vals[beta_idx] - gamma_ref[closest_beta_ref_idx])
                    omega_diff.append(omega_vals[beta_idx] - omega_ref[closest_beta_ref_idx])
            ax1.plot(beta_to_compare, omega_diff, label=labels[folder_idx], marker=marker_list[folder_idx], mfc="none", markersize=marker_size)
            ax2.plot(beta_to_compare, gamma_diff, label=labels[folder_idx], marker=marker_list[folder_idx], mfc="none", markersize=marker_size)


    for ax in [ax1, ax2]:
        ax.grid(True)

    ax1.legend(loc="best")
    if not omega_diff:
        ax1.set_ylabel(r"$\omega$")
        ax2.set_ylabel(r"$\gamma$")
    else:
        ax1.set_ylabel(r"$\omega-\omega_{ref}$")
        ax2.set_ylabel(r"$\gamma-\gamma_{ref}$")
    ax2.set_xlabel(r"$\beta$")
    plt.show()

    return

def make_beta_pair_plots_from_pickles(pickle_longname_pairs, labels, marker_size=1, omega_diff=False):
    """Either plot Omega(beta), or, if omega_diff=True, plot (Omega-Omega_ref) (beta),
    where Omega_ref is taken from the first pickle. """

    marker_list = ["s", "o", "P", "X", "v", "^", "<", ">", "1", "2", "3"]
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)

    for pair_idx, pickle_longname_pair in enumerate(pickle_longname_pairs):
        pickle_longname_1 = pickle_longname_pair[0]
        pickle_longname_2 = pickle_longname_pair[1]
        myfile = open(pickle_longname_1, "rb")
        [beta_vals, gamma_vals, omega_vals] = pickle.load(myfile)
        myfile.close()
        # This is the reference beta
        beta_ref = beta_vals
        gamma_ref = gamma_vals
        omega_ref = omega_vals

        myfile = open(pickle_longname_2, "rb")
        [beta_vals, gamma_vals, omega_vals] = pickle.load(myfile)
        myfile.close()
        # Find where beta vals match.
        beta_to_compare = []
        gamma_diff = []
        omega_diff = []
        for beta_idx, beta_val in enumerate(beta_vals):
            print("beta_val, beta_ref = ", beta_val, beta_ref)
            ## Find if the beta val is within a small tolerance of any
            ## beta_ref vals
            closest_beta_ref_idx = np.argmin(abs(beta_ref - beta_val))
            if abs(beta_ref[closest_beta_ref_idx] - beta_val) < 1e-6:
                beta_to_compare.append(beta_val)
                gamma_diff.append(gamma_vals[beta_idx] - gamma_ref[closest_beta_ref_idx])
                omega_diff.append(omega_vals[beta_idx] - omega_ref[closest_beta_ref_idx])
        ax1.plot(beta_to_compare, omega_diff, label=labels[pair_idx], marker=marker_list[pair_idx], mfc="none", markersize=marker_size)
        ax2.plot(beta_to_compare, gamma_diff, label=labels[pair_idx], marker=marker_list[pair_idx], mfc="none", markersize=marker_size)


    for ax in [ax1, ax2]:
        ax.grid(True)

    ax1.legend(loc="best")
    if not omega_diff:
        ax1.set_ylabel(r"$\omega$")
        ax2.set_ylabel(r"$\gamma$")
    else:
        ax1.set_ylabel(r"$\omega-\omega_{ref}$")
        ax2.set_ylabel(r"$\gamma-\gamma_{ref}$")
    ax2.set_xlabel(r"$\beta$")
    plt.show()

    return

def plot_beta_scans_with_resolution_checks():
    """ """
    ## ky=0.5. Trying a few things out
    make_beta_plots_from_pickles(
                [
                "sims/" + gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + gs2_beta_scan_ky_05_np3_nt48_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np2_nt64_nvpa24_nmu18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np2_nt64_nvpa24_nmu18_fapar1_fbpar1_streaming_implicit_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_streaming_implicit_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_streaming_implicit_folder + "/beta_gamma_omega.pickle",
                ],
                [
                "GS2, hr",
                "GS2, lr",
                "stella, hr",
                "stella, lr + nt64",
                "stella, vlr + nt64",
                "stella, vlr + nt32",
                "stella, lr (str. impl.) + nt64",
                "stella, vlr (str. impl.) + nt64",
                "stella, vlr (str. impl.) + nt32",
                ],
                marker_size=8,
                omega_diff=True
                )
    ### Benchmark stella with implicit streaming
    make_beta_plots_from_pickles(
                [
                "sims/" + gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + gs2_beta_scan_ky_05_np3_nt48_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np2_nt64_nvpa24_nmu18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/viking_gfortran_build/" + stella_beta_scan_ky_05_np2_nt64_nvpa24_nmu18_fapar1_fbpar1_streaming_implicit_folder + "/beta_gamma_omega.pickle",
                "sims/viking_gfortran_build/" + stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_streaming_implicit_folder + "/beta_gamma_omega.pickle",
                "sims/viking_gfortran_build/" + stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_streaming_implicit_folder + "/beta_gamma_omega.pickle",
                ],
                [
                "GS2, hr",
                "GS2, lr",
                "stella, hr",
                "stella, lr + nt64",
                "stella, vlr + nt64",
                "stella, vlr + nt32",
                "stella, lr (str. impl.) + nt64",
                "stella, vlr (str. impl.) + nt64",
                "stella, vlr (str. impl.) + nt32",
                ],
                marker_size=8,
                omega_diff=True
                )

    # make_beta_plots_from_pickles(
    #             [
    #             "sims/" + gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + gs2_beta_scan_ky_05_np3_nt48_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + stella_beta_scan_ky_05_np2_nt64_nvpa24_nmu18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             ],
    #             [
    #             "GS2, hr",
    #             "GS2, lr",
    #             "stella, hr",
    #             "stella, lr + nt64",
    #             "stella, vlr + nt64",
    #             "stella, vlr + nt32",
    #             ],
    #             marker_size=8,
    #             omega_diff=True
    #             )

    # make_beta_plots_from_pickles(
    #             [
    #             "sims/" + old_stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + new_stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + old_stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
    #              "sims/" + new_stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + old_stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
    #              "sims/" + new_stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
    #             ],
    #             [
    #             "electromagnetic vlr + nt128",
    #             "electromagnetic-new vlr + nt128",
    #             "electromagnetic vlr + nt256",
    #             "electromagnetic-new vlr + nt256",
    #             "electromagnetic vlr + nt512",
    #             "electromagnetic-new vlr + nt512",
    #             ],
    #             marker_size=8,
    #             )

    # make_beta_pair_plots_from_pickles(
    #             [
    #             ["sims/" + old_stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
    #              "sims/" + new_stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle"],
    #             ["sims/" + old_stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
    #              "sims/" + new_stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle"],
    #             ["sims/" + old_stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
    #              "sims/" + new_stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle"],
    #             ],
    #             [
    #             "vlr + nt128",
    #             "vlr + nt256",
    #             "vlr + nt512",
    #             ],
    #             marker_size=8,
    #             )

    ## ky=0.5. ntheta resolution
    # make_beta_plots_from_pickles(
    #             [
    #             "sims/" + gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + gs2_beta_scan_ky_05_np3_nt48_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
    #             ],
    #             [
    #             "GS2, hr",
    #             "GS2, lr",
    #             "stella, vlr + nt64",
    #             "stella, vlr + nt128",
    #             "stella, vlr + nt256",
    #             "stella, vlr + nt512",
    #             ],
    #             marker_size=8,
    #             omega_diff=True
    #             )


    ## ky=1. ntheta resolution
    # make_beta_plots_from_pickles(
    #             [
    #             "sims/" + gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_ky_1_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + gs2_beta_scan_ky_05_np3_nt48_ng8_ne18_ky_1_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_ky_1_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_ky_1_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_ky_1_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_ky_1_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             ],
    #             [
    #             "GS2, hr",
    #             "GS2, lr",
    #             "stella, vlr + nt64",
    #             "stella, vlr + nt128",
    #             "stella, vlr + nt256",
    #             "stella, vlr + nt512",
    #             ],
    #             marker_size=8
    #             )

    ## ky=5. ntheta resolution
    # make_beta_plots_from_pickles(
    #             [
    #             "sims/" + gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_ky_5_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + gs2_beta_scan_ky_05_np3_nt48_ng8_ne18_ky_5_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_ky_5_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_ky_5_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_ky_5_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             "sims/" + stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_ky_5_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
    #             ],
    #             [
    #             "GS2, hr",
    #             "GS2, lr",
    #             "stella, vlr + nt64",
    #             "stella, vlr + nt128",
    #             "stella, vlr + nt256",
    #             "stella, vlr + nt512",
    #             ],
    #             marker_size=8
    #             )

    return

def make_beta_scan_plots_for_ipp_greifswald_talk():
    """ """

    def make_beta_plots_from_pickles_for_talk(pickle_longnames, labels,
                                        marker_size=1, lw=1, omega_diff=False):
        """Either plot Omega(beta), or, if omega_diff=True, plot (Omega-Omega_ref) (beta),
        where Omega_ref is taken from the first pickle. """

        marker_list = ["s", "o", "P", "X", "v", "^", "<", ">", "1", "2", "3"]
        fig = plt.figure(figsize=(12,8))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex=ax1)

        for folder_idx, pickle_longname in enumerate(pickle_longnames):
            myfile = open(pickle_longname, "rb")
            [beta_vals, gamma_vals, omega_vals] = pickle.load(myfile)
            myfile.close()
            if not omega_diff:
                ax1.plot(beta_vals, omega_vals, label=labels[folder_idx], marker=marker_list[folder_idx], markersize=marker_size, lw=lw, mfc="none",)
                ax2.plot(beta_vals, gamma_vals, label=labels[folder_idx], marker=marker_list[folder_idx], markersize=marker_size, lw=lw, mfc="none",)
            else:
                if folder_idx==0:
                    # This is the reference beta
                    beta_ref = beta_vals
                    gamma_ref = gamma_vals
                    omega_ref = omega_vals
                ## Whether or not this is the reference beta scan, plot Omega-Omega_ref
                # Find where beta vals match.
                beta_to_compare = []
                gamma_diff = []
                omega_diff = []
                for beta_idx, beta_val in enumerate(beta_vals):
                    ## Find if the beta val is within a small tolerance of any
                    ## beta_ref vals
                    closest_beta_ref_idx = np.argmin(abs(beta_ref - beta_val))
                    if abs(beta_ref[closest_beta_ref_idx] - beta_val) < 1e-6:
                        beta_to_compare.append(beta_val)
                        gamma_diff.append(gamma_vals[beta_idx] - gamma_ref[closest_beta_ref_idx])
                        omega_diff.append(omega_vals[beta_idx] - omega_ref[closest_beta_ref_idx])
                ax1.plot(beta_to_compare, np.array(omega_diff)*100, label=labels[folder_idx], marker=marker_list[folder_idx], mfc="none", markersize=marker_size, lw=lw)
                ax2.plot(beta_to_compare, np.array(gamma_diff)*100, label=labels[folder_idx], marker=marker_list[folder_idx], mfc="none", markersize=marker_size, lw=lw)

        label_fontsize = 14

        for ax in [ax1, ax2]:
            ax.grid(True)
            ax.tick_params("x", which="both", labelsize=label_fontsize)
            ax.tick_params("y", which="both", labelsize=label_fontsize)
            #ax.set_yticks(fontsize=label_fontsize)
        ax1.tick_params("x", which="both", labelbottom=False)
        my_legendfontsize = 14
        my_fontsize = 20
        my_rotation=70
        my_padding = 20
        ax1.legend(loc="best", fontsize=my_legendfontsize)
        if not omega_diff:
            ax1.set_ylabel(r"frequency $\omega$", fontsize=my_fontsize, rotation=my_rotation, labelpad=my_padding)
            ax2.set_ylabel(r"growth rate $\gamma$", fontsize=my_fontsize, rotation=my_rotation, labelpad=my_padding)
            save_name = "for_ipp_a.png"
        else:
            save_name = "for_ipp_b.png"
            ax1.set_ylabel(r"$\Delta\omega (\%)$", fontsize=my_fontsize, labelpad=my_padding, rotation=my_rotation)
            ax2.set_ylabel(r"$\Delta\gamma (\%)$", fontsize=my_fontsize, labelpad=my_padding, rotation=my_rotation)
        ax2.set_xlabel(r"$\beta$", fontsize=my_fontsize)
        plt.tight_layout()
        plt.show()
        # plt.savefig(save_name)
        # plt.close()

        return
    ## ky=0.5. Trying a few things out
    make_beta_plots_from_pickles_for_talk(
                [
                "sims/" + gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_bespoke_explicit_folder + "/beta_gamma_omega.pickle",
                # "sims/" + stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_bespoke_stream_implicit_folder + "/beta_gamma_omega.pickle",
                # "sims/" + stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_streaming_implicit_folder + "/beta_gamma_omega.pickle",
                ],
                [
                "GS2",
                "stella, fully explicit (RK3)",
                # "stella, vlr + nt32",
                "stella, streaming implicit",
                # "stella, vlr (str. impl.) + nt32",
                ],
                marker_size=10,
                lw=3,
                omega_diff=False
                )
    make_beta_plots_from_pickles_for_talk(
                [
                "sims/" + gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_bespoke_explicit_folder + "/beta_gamma_omega.pickle",
                # "sims/" + stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_bespoke_stream_implicit_folder + "/beta_gamma_omega.pickle",
                # "sims/" + stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_streaming_implicit_folder + "/beta_gamma_omega.pickle",
                ],
                [
                "GS2",
                "stella, fully explicit (RK3)",
                # "stella, vlr + nt32",
                "stella, streaming implicit",
                # "stella, vlr (str. impl.) + nt32",
                ],
                marker_size=10,
                lw=3,
                omega_diff=True
                )
    ### Benchmark stella with implicit streaming

    return

def plot_beta_scans_with_resolution_checks_for_tdotp():
    """ """

    def make_beta_plots_from_pickles_local(pickle_longnames, labels, save_name,
                                        marker_size=1, lw=1, omega_diff=False,
                                        fractional=False,
                                        omega_diff_ylabel=r"$\omega - \omega_{ref}$",
                                        gamma_diff_ylabel=r"$\gamma - \gamma_{ref}$",
                                        ylabel_fontsize=14,
                                        xlabel_fontsize=14
                                        ):
        """Either plot Omega(beta), or, if omega_diff=True, plot (Omega-Omega_ref) (beta),
        where Omega_ref is taken from the first pickle. """

        marker_list = ["s", "o", "P", "X", "v", "^", "<", ">", "1", "2", "3"]
        fig = plt.figure(figsize=(12,8))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex=ax1)

        for folder_idx, pickle_longname in enumerate(pickle_longnames):
            myfile = open(pickle_longname, "rb")
            [beta_vals, gamma_vals, omega_vals] = pickle.load(myfile)
            myfile.close()
            if not omega_diff:
                ax1.plot(beta_vals, omega_vals, label=labels[folder_idx], marker=marker_list[folder_idx], mfc="none", markersize=marker_size)
                ax2.plot(beta_vals, gamma_vals, label=labels[folder_idx], marker=marker_list[folder_idx], mfc="none", markersize=marker_size)
            else:
                if folder_idx==0:
                    # This is the reference beta
                    beta_ref = beta_vals
                    gamma_ref = gamma_vals
                    omega_ref = omega_vals
                ## Whether or not this is the reference beta scan, plot Omega-Omega_ref
                # Find where beta vals match.
                beta_to_compare = []
                gamma_diff = []
                omega_diff = []
                for beta_idx, beta_val in enumerate(beta_vals):
                    ## Find if the beta val is within a small tolerance of any
                    ## beta_ref vals
                    closest_beta_ref_idx = np.argmin(abs(beta_ref - beta_val))
                    if abs(beta_ref[closest_beta_ref_idx] - beta_val) < 1e-6:
                        beta_to_compare.append(beta_val)
                        if fractional:
                            gamma_diff.append((gamma_vals[beta_idx] - gamma_ref[closest_beta_ref_idx])*100/gamma_ref[closest_beta_ref_idx])
                            omega_diff.append((omega_vals[beta_idx] - omega_ref[closest_beta_ref_idx])*100/omega_ref[closest_beta_ref_idx])
                        else:
                            gamma_diff.append(gamma_vals[beta_idx] - gamma_ref[closest_beta_ref_idx])
                            omega_diff.append(omega_vals[beta_idx] - omega_ref[closest_beta_ref_idx])
                ax1.plot(beta_to_compare, omega_diff, label=labels[folder_idx], marker=marker_list[folder_idx], mfc="none", markersize=marker_size)
                ax2.plot(beta_to_compare, gamma_diff, label=labels[folder_idx], marker=marker_list[folder_idx], mfc="none", markersize=marker_size)


        for ax in [ax1, ax2]:
            ax.grid(True)

        ax1.legend(loc="best")
        if not omega_diff:
            ax1.set_ylabel(r"$\omega$", fontsize=ylabel_fontsize)
            ax2.set_ylabel(r"$\gamma$", fontsize=ylabel_fontsize)
        else:
            ax1.set_ylabel(omega_diff_ylabel, fontsize=ylabel_fontsize)
            ax2.set_ylabel(gamma_diff_ylabel, fontsize=ylabel_fontsize)
        ax2.set_xlabel(r"$\beta$", fontsize=xlabel_fontsize)
        plt.savefig(save_name)
        plt.close()

        return

    ## ky=0.5. Trying a few things out
    make_beta_plots_from_pickles_local(
                [
                "sims/" + gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + gs2_beta_scan_ky_05_np3_nt48_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + new_stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
                "sims/" + new_stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
                "sims/" + new_stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2_streaming_implicit_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2_streaming_implicit_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2_streaming_implicit_folder + "/beta_gamma_omega.pickle",
                ],
                [
                "GS2, hr (ntheta=64)",
                "GS2, lr (ntheta=48)",
                "stella, hr (ntheta=64)",
                "stella, lr (ntheta=64)",
                "stella, lr (ntheta=128)",
                "stella, lr (ntheta=256)",
                "stella, lr (ntheta=512)",
                "stella, lr, str impl (ntheta=128)",
                "stella, lr, str impl (ntheta=256)",
                "stella, lr, str impl (ntheta=512)",
                ],
                "absolute_comparison.png",
                marker_size=8, lw=2,
                ylabel_fontsize=14,
                xlabel_fontsize=14
                )

    make_beta_plots_from_pickles_local(
                [
                "sims/" + gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + gs2_beta_scan_ky_05_np3_nt48_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + new_stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
                "sims/" + new_stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
                "sims/" + new_stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
                ],
                [
                "GS2, hr",
                "GS2, lr",
                "stella, hr",
                "stella, lr + nt64",
                "stella, lr + nt128",
                "stella, lr + nt256",
                "stella, lr + nt512",
                ],
                "relative_comparison.png",
                marker_size=8, lw=2,
                omega_diff=True,
                omega_diff_ylabel=r"$\omega - \omega_{GS2, hr}$",
                gamma_diff_ylabel=r"$\gamma - \gamma_{GS2, hr}$",
                ylabel_fontsize=14,
                xlabel_fontsize=14
                )

    make_beta_plots_from_pickles_local(
                [
                "sims/" + gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + gs2_beta_scan_ky_05_np3_nt48_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + new_stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
                "sims/" + new_stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
                "sims/" + new_stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2_folder + "/beta_gamma_omega.pickle",
                ],
                [
                "GS2, hr (ntheta=64)",
                "GS2, lr (ntheta=48)",
                "stella, hr (ntheta=64)",
                "stella, lr (ntheta=64)",
                "stella, lr (ntheta=128)",
                "stella, lr (ntheta=256)",
                "stella, lr (ntheta=512)",
                ],
                "fractional_relative_comparison.png",
                marker_size=8, lw=2,
                omega_diff=True,
                omega_diff_ylabel=r"$(\omega - \omega_{GS2, hr})/\omega_{GS2, hr}$ (%)",
                gamma_diff_ylabel=r"$(\gamma - \gamma_{GS2, hr})/\omega_{GS2, hr}$ (%)",
                ylabel_fontsize=14,
                xlabel_fontsize=14,
                fractional=True
                )


    return

if __name__ == "__main__":
    # plot_different_beta_scans()
    # plot_stella_scan_vs_gs2_pickle()
    # plot_stella_scan_vs_gs2_pickle_for_poster()
    # plot_beta_scans_with_resolution_checks()
    # plot_beta_scans_with_resolution_checks_for_tdotp()
    make_beta_scan_plots_for_ipp_greifswald_talk()
