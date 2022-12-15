""" """

import sys
sys.path.append("../postprocessing_tools")
from helper_ncdf_new import view_ncdf_variables, extract_data_from_ncdf
import matplotlib.pyplot as plt
import numpy as np
import glob
import re
import pickle
import make_param_scans as mps

default_cmap = plt.get_cmap("tab10")

stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1"
stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2"
stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2"
stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2"
stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1_folder="stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1"

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

def make_beta_plots_from_pickles_for_thesis_src_h(pickle_longnames, labels, marker_size=1, omega_diff=False,
                                save_name="test.eps"):
    """Either plot Omega(beta), or, if omega_diff=True, plot (Omega-Omega_ref) (beta),
    where Omega_ref is taken from the first pickle. """

    marker_size = 12
    label_fontsize = 40
    legend_fontsize = 14
    marker_list = ["s", "o", "P", "X", "v", "^", "<", ">", "1", "2", "3"]
    lw_list = [4, 3, 3, 2,
                1, 1, 1, 1, 1, 1, 1,]
    ls_list = ["-", "--", "-.", (0, (4,1,2,1)) ,
                "-", "--", "-.", (0, (4,1,2,1))
                ]

    my_xticklength = 7
    my_xtickwidth = 2
    my_xticklength_minor = 4
    my_xtickwidth_minor = 1

    x_ticklabelfontsize = 20
    y_ticklabelfontsize = 20
    y_labelfontsize = 30
    x_labelfontsize = 30
    bracket_fontsize = 70
    itg_fontsize = 30

    top = 0.98
    left = 0.14
    right = 0.98
    bottom = 0.13
    vspace = 0.02
    height = (top - bottom - vspace)/2
    row1_bottom = bottom + height + vspace
    width = right - left

    fig = plt.figure(figsize=(12,12))
    ax1 = fig.add_axes((left, row1_bottom, width, height))
    ax2 = fig.add_axes((left, bottom, width, height))

    ## Use the second folder as the "correct" one
    reference_pickle_longname = pickle_longnames[1]
    myfile = open(reference_pickle_longname, "rb")
    [beta_ref, gamma_ref, omega_ref] = pickle.load(myfile)
    myfile.close()
    col_list = [default_cmap(0), default_cmap(1), default_cmap(2), default_cmap(3),
    default_cmap(4), default_cmap(5), default_cmap(6), default_cmap(7),
    default_cmap(8), default_cmap(9), default_cmap(10), default_cmap(11),
            ]

    for folder_idx, pickle_longname in enumerate(pickle_longnames):
        myfile = open(pickle_longname, "rb")
        [beta_vals, gamma_vals, omega_vals] = pickle.load(myfile)
        myfile.close()
        if not omega_diff:
            ax1.plot(beta_vals, omega_vals, label=labels[folder_idx],
                lw=lw_list[folder_idx], ls=ls_list[folder_idx],
                c=col_list[folder_idx],
                marker=marker_list[folder_idx], mfc="none", markersize=marker_size)
            ax2.plot(beta_vals, gamma_vals, label=labels[folder_idx],
                lw=lw_list[folder_idx], ls=ls_list[folder_idx],
                c=col_list[folder_idx],
                marker=marker_list[folder_idx], mfc="none", markersize=marker_size)
        else:
            if folder_idx!=1:
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
                        gamma_diff.append(100*(gamma_vals[beta_idx] - gamma_ref[closest_beta_ref_idx])/gamma_ref[closest_beta_ref_idx])
                        omega_diff.append(100*(omega_vals[beta_idx] - omega_ref[closest_beta_ref_idx])/omega_ref[closest_beta_ref_idx])
                ax1.plot(beta_to_compare, omega_diff, label=labels[folder_idx],
                         lw=lw_list[folder_idx], ls=ls_list[folder_idx],
                         c=col_list[folder_idx],
                         marker=marker_list[folder_idx], mfc="none", markersize=marker_size)
                ax2.plot(beta_to_compare, gamma_diff, label=labels[folder_idx],
                         lw=lw_list[folder_idx], ls=ls_list[folder_idx],
                         c=col_list[folder_idx],
                         marker=marker_list[folder_idx], mfc="none", markersize=marker_size)

    # for ax in [ax1, ax2]:
    #     ax.grid(True)


    ax1.legend(loc="best", fontsize=legend_fontsize)
    min_beta = -0.001
    max_beta = 0.041
    for ax in [ax1, ax2]:
        ax.set_xlim((min_beta, max_beta))
        ax.tick_params("y", which="major", length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.set_xticks([0, 0.01, 0.02, 0.03, 0.04])
        ax.set_xticklabels([r"$0$", r"$0.01$", r"$0.02$", r"$0.03$", r"$0.04$"], fontsize=x_ticklabelfontsize)
    ax1.tick_params("x", which="major", length=my_xticklength, width=my_xtickwidth, direction="in", labelbottom=False)
    ax2.tick_params("x", which="major", length=my_xticklength, width=my_xtickwidth, direction="in", labelbottom=True)

    if not omega_diff:
        min_freq = 0.16
        max_freq = 0.92
        min_gamma = 0.01
        max_gamma = 0.79
        ax1.set_ylim((min_freq, max_freq))
        ax2.set_ylim((min_gamma, max_gamma))
        ax1.set_yticks([0.2, 0.4, 0.6, 0.8])
        ax1.set_yticklabels([r"$0.2$", r"$0.4$", r"$0.6$", r"$0.8$"], fontsize=x_ticklabelfontsize)
        ax2.set_yticks([0.2, 0.4, 0.6])
        ax2.set_yticklabels([r"$0.2$", r"$0.4$", r"$0.6$"], fontsize=x_ticklabelfontsize)
        ax1.set_ylabel(r"$\tilde{\omega}$", fontsize=label_fontsize)
        ax2.set_ylabel(r"$\tilde{\gamma}$", fontsize=label_fontsize)
        # ax1.annotate(r"$\{$", (0.002, 0.4), rotation=90, fontsize=bracket_fontsize)
        itg_x = 0.005
        kbm_x = 0.027
        ax1.annotate(r"ITG", xy=(itg_x, 0.4), xytext=(itg_x, 0.43), fontsize=itg_fontsize,
                    ha='center', va='bottom',
                        arrowprops=dict(arrowstyle='-[, widthB=3.0, lengthB=0.5', lw=2.0) )
        ax1.annotate(r"KBM", xy=(kbm_x, 0.38), xytext=(kbm_x, 0.35), fontsize=itg_fontsize,
                    ha='center', va='top',
                        arrowprops=dict(arrowstyle='-[, widthB=7.0, lengthB=0.5', lw=2.0) )


    else:
        ax1.plot((min_beta, max_beta), (0,0), c="gray", zorder=-10)
        ax2.plot((min_beta, max_beta), (0,0), c="gray", zorder=-10)
        ax1.set_ylabel(r"$\Delta\tilde{\omega}_{ref}$ ($\%$)", fontsize=label_fontsize)
        ax2.set_ylabel(r"$\Delta\tilde{\gamma}_{ref}$ ($\%$)", fontsize=label_fontsize)
        ax1.set_yticks([0, 2, 4, 6, 8])
        ax1.set_yticklabels([r"$0$", r"$2$", r"$4$", r"$6$", r"$8$"], fontsize=x_ticklabelfontsize)
        ax2.set_yticks([-20, -10, 0, 10])
        ax2.set_yticklabels([r"$-20$", r"$-10$", r"$0$", r"$10$"], fontsize=x_ticklabelfontsize)
    ax2.set_xlabel(r"$\beta$", fontsize=label_fontsize)
    plt.savefig(save_name)
    plt.close()

    return

def make_beta_plots_from_pickles_for_thesis_gbar(pickle_longnames, labels, marker_size=1, omega_diff=False,
                                save_name="test.eps", for_talk=False):
    """Either plot Omega(beta), or, if omega_diff=True, plot (Omega-Omega_ref) (beta),
    where Omega_ref is taken from the first pickle. """

    marker_size = 12
    label_fontsize = 40
    legend_fontsize = 14
    marker_list = ["s", "o", "P", "X", "v", "^", "<", ">", "1", "2", "3"]
    lw_list = [4, 3, 3, 2,
                2, 2, 2, 2, 1, 1, 1,]
    ls_list = ["-", "--", "-.", (0, (4,1,2,1)) ,
                "-", "--", "-.", (0, (4,1,2,1))
                ]

    my_xticklength = 7
    my_xtickwidth = 2
    my_xticklength_minor = 4
    my_xtickwidth_minor = 1

    x_ticklabelfontsize = 20
    y_ticklabelfontsize = 20
    y_labelfontsize = 30
    x_labelfontsize = 30
    bracket_fontsize = 70
    itg_fontsize = 30

    top = 0.98
    left = 0.14
    right = 0.98
    bottom = 0.13
    vspace = 0.02
    height = (top - bottom - vspace)/2
    row1_bottom = bottom + height + vspace
    width = right - left

    fig = plt.figure(figsize=(12,12))
    ax1 = fig.add_axes((left, row1_bottom, width, height))
    ax2 = fig.add_axes((left, bottom, width, height))


    ## Use the second folder as the "correct" one
    reference_pickle_longname = pickle_longnames[1]
    myfile = open(reference_pickle_longname, "rb")
    [beta_ref, gamma_ref, omega_ref] = pickle.load(myfile)
    myfile.close()
    col_list = [default_cmap(0), default_cmap(1), default_cmap(2), default_cmap(3),
    default_cmap(4), default_cmap(5), default_cmap(6), default_cmap(7),
    default_cmap(8), default_cmap(9), default_cmap(10), default_cmap(11),
            ]

    for folder_idx, pickle_longname in enumerate(pickle_longnames):
        myfile = open(pickle_longname, "rb")
        [beta_vals, gamma_vals, omega_vals] = pickle.load(myfile)
        myfile.close()
        if not omega_diff:
            ax1.plot(beta_vals, omega_vals, label=labels[folder_idx],
                lw=lw_list[folder_idx], ls=ls_list[folder_idx],
                c=col_list[folder_idx],
                marker=marker_list[folder_idx], mfc="none", markersize=marker_size)
            ax2.plot(beta_vals, gamma_vals, label=labels[folder_idx],
                lw=lw_list[folder_idx], ls=ls_list[folder_idx],
                c=col_list[folder_idx],
                marker=marker_list[folder_idx], mfc="none", markersize=marker_size)
        else:
            if folder_idx!=1:
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
                        gamma_diff.append(100*(gamma_vals[beta_idx] - gamma_ref[closest_beta_ref_idx])/gamma_ref[closest_beta_ref_idx])
                        omega_diff.append(100*(omega_vals[beta_idx] - omega_ref[closest_beta_ref_idx])/omega_ref[closest_beta_ref_idx])
                print("sim:", labels[folder_idx])
                print("beta_to_compare = ", beta_to_compare)
                print("gamma_diff = ", gamma_diff)
                print("omega_diff = ", omega_diff)
                ax1.plot(beta_to_compare, omega_diff, label=labels[folder_idx],
                         lw=lw_list[folder_idx], ls=ls_list[folder_idx],
                         c=col_list[folder_idx],
                         marker=marker_list[folder_idx], mfc="none", markersize=marker_size)
                ax2.plot(beta_to_compare, gamma_diff, label=labels[folder_idx],
                         lw=lw_list[folder_idx], ls=ls_list[folder_idx],
                         c=col_list[folder_idx],
                         marker=marker_list[folder_idx], mfc="none", markersize=marker_size)

    # for ax in [ax1, ax2]:
    #     ax.grid(True)

    if for_talk:
        ax1.legend(loc="upper right", fontsize=legend_fontsize)
    else:
        if omega_diff:
            ax1.legend(loc="middle right", fontsize=legend_fontsize)
        else:
            ax1.legend(loc="upper right", fontsize=legend_fontsize)
    min_beta = -0.001
    max_beta = 0.041
    for ax in [ax1, ax2]:
        ax.set_xlim((min_beta, max_beta))
        ax.tick_params("y", which="major", length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.set_xticks([0, 0.01, 0.02, 0.03, 0.04])
        ax.set_xticklabels([r"$0$", r"$0.01$", r"$0.02$", r"$0.03$", r"$0.04$"], fontsize=x_ticklabelfontsize)
    ax1.tick_params("x", which="major", length=my_xticklength, width=my_xtickwidth, direction="in", labelbottom=False)
    ax1.tick_params("y", which="minor", length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in", labelbottom=False)
    ax2.tick_params("x", which="major", length=my_xticklength, width=my_xtickwidth, direction="in", labelbottom=True)
    ax2.tick_params("y", which="minor", length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in", labelbottom=False)

    if not omega_diff:
        min_freq = 0.16
        max_freq = 1.5
        min_gamma = 0.01
        max_gamma = 0.79
        ax1.set_ylim((min_freq, max_freq))
        ax2.set_ylim((min_gamma, max_gamma))
        ax1.set_yticks([0.4, 0.8, 1.2])
        ax1.set_yticks([0.2, 0.6, 1, 1.4], minor=True)
        ax1.set_yticklabels([r"$0.4$", r"$0.8$", r"$1.2$"], fontsize=x_ticklabelfontsize)
        ax2.set_yticks([0.2, 0.4, 0.6])
        ax2.set_yticklabels([r"$0.2$", r"$0.4$", r"$0.6$"], fontsize=x_ticklabelfontsize)
        ax1.set_ylabel(r"$\tilde{\omega}$", fontsize=label_fontsize)
        ax2.set_ylabel(r"$\tilde{\gamma}$", fontsize=label_fontsize)
        # ax1.annotate(r"$\{$", (0.002, 0.4), rotation=90, fontsize=bracket_fontsize)
        itg_x = 0.005
        kbm_x = 0.027
        ax1.annotate(r"ITG", xy=(itg_x, 0.4), xytext=(itg_x, 0.43), fontsize=itg_fontsize,
                    ha='center', va='bottom',
                        arrowprops=dict(arrowstyle='-[, widthB=3.0, lengthB=0.5', lw=2.0) )
        ax1.annotate(r"KBM", xy=(kbm_x, 0.38), xytext=(kbm_x, 0.35), fontsize=itg_fontsize,
                    ha='center', va='top',
                        arrowprops=dict(arrowstyle='-[, widthB=7.0, lengthB=0.5', lw=2.0) )

    else:
        min_freq = -25
        max_freq = 20
        min_gamma = -100
        max_gamma = 45
        ax1.set_ylim((min_freq, max_freq))
        ax2.set_ylim((min_gamma, max_gamma))
        ax1.plot((min_beta, max_beta), (0,0), c="gray", zorder=-10)
        ax2.plot((min_beta, max_beta), (0,0), c="gray", zorder=-10)
        ax1.set_ylabel(r"$\Delta\tilde{\omega}_{ref}$ ($\%$)", fontsize=label_fontsize)
        ax2.set_ylabel(r"$\Delta\tilde{\gamma}_{ref}$ ($\%$)", fontsize=label_fontsize)
        ax1.set_yticks([-20, -10, 0, 10, 20])
        ax1.set_yticks([-15, -5, 5, 15], minor=True)
        ax1.set_yticklabels([r"$-20$", r"$-10$", r"$0$", r"$10$", r"$20$"], fontsize=x_ticklabelfontsize)
        ax2.set_yticks([-100, -50, 0])
        ax2.set_yticks([-75, -25, 25], minor=True)
        ax2.set_yticklabels([r"$-100$", r"$-50$", r"$0$"], fontsize=x_ticklabelfontsize)
    ax2.set_xlabel(r"$\beta$", fontsize=label_fontsize)
    plt.savefig(save_name)
    plt.close()

    return

def make_omega_beta_plots_fapar1_fbpar1():
    """Plot CBC beta scan for thesis, showing omega and gamma,
    rather than the difference between these """

    make_beta_plots_from_pickles_for_thesis_src_h(
                [
                "sims/" + mps.gs2_beta_scan_ky_05_np2_nt32_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + mps.gs2_beta_scan_ky_05_np4_nt64_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + mps.stella_src_h_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_zupw0_fapar1_fbpar1_str_mirror_implicit_folder + "/beta_gamma_omega.pickle",
                "sims/" + mps.stella_src_h_beta_scan_ky_05_np4_nt64_nvpa18_nmu12_zupw0_fapar1_fbpar1_str_mirror_implicit_folder + "/beta_gamma_omega.pickle",
                ],
                [
                r"GS2 ($n_z=32$, $n_{period}=2$)",
                r"GS2 ($n_z=64$, $n_{period}=4$)",
                r"stella (implicit) ($n_z=32$, $n_{period}=2$)",
                r"stella (implicit) ($n_z=64$, $n_{period}=4$)",
                ],
                omega_diff=False,
                save_name="cbc_Omega_beta_scan.eps"

                                )
    make_beta_plots_from_pickles_for_thesis_src_h(
                [
                "sims/" + mps.gs2_beta_scan_ky_05_np2_nt32_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + mps.gs2_beta_scan_ky_05_np4_nt64_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + mps.stella_src_h_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_zupw0_fapar1_fbpar1_str_mirror_implicit_folder + "/beta_gamma_omega.pickle",
                "sims/" + mps.stella_src_h_beta_scan_ky_05_np4_nt64_nvpa18_nmu12_zupw0_fapar1_fbpar1_str_mirror_implicit_folder + "/beta_gamma_omega.pickle",
                ],
                [
                r"GS2 ($n_z=32$, $n_{period}=2$)",
                r"Gs2 ($n_z=64$, $n_{period}=4$)",
                r"stella (implicit) ($n_z=32$, $n_{period}=2$)",
                r"stella (implicit) ($n_z=64$, $n_{period}=4$)",
                ],
                omega_diff=True,
                save_name="cbc_Omegadiff_beta_scan.eps"

                                )
    return

def for_talk_make_omega_beta_plots_fapar1_fbpar1():
    """Plot CBC beta scan for thesis, showing omega and gamma,
    rather than the difference between these """

    make_beta_plots_from_pickles_for_thesis_src_h(
                [
                "sims/" + mps.gs2_beta_scan_ky_05_np2_nt32_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + mps.gs2_beta_scan_ky_05_np4_nt64_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + mps.stella_src_h_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_zupw0_fapar1_fbpar1_str_mirror_implicit_folder + "/beta_gamma_omega.pickle",
                "sims/" + mps.stella_src_h_beta_scan_ky_05_np4_nt64_nvpa18_nmu12_zupw0_fapar1_fbpar1_str_mirror_implicit_folder + "/beta_gamma_omega.pickle",
                ],
                [
                r"GS2 ($n_z=32$, $n_{period}=2$)",
                r"GS2 ($n_z=64$, $n_{period}=4$)",
                r"stella ($n_z=32$, $n_{period}=2$)",
                r"stella ($n_z=64$, $n_{period}=4$)",
                ],
                omega_diff=False,
                save_name="for_talk_cbc_Omega_beta_scan.eps"

                                )
    make_beta_plots_from_pickles_for_thesis_src_h(
                [
                "sims/" + mps.gs2_beta_scan_ky_05_np2_nt32_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + mps.gs2_beta_scan_ky_05_np4_nt64_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + mps.stella_src_h_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_zupw0_fapar1_fbpar1_str_mirror_implicit_folder + "/beta_gamma_omega.pickle",
                "sims/" + mps.stella_src_h_beta_scan_ky_05_np4_nt64_nvpa18_nmu12_zupw0_fapar1_fbpar1_str_mirror_implicit_folder + "/beta_gamma_omega.pickle",
                ],
                [
                r"GS2 ($n_z=32$, $n_{period}=2$)",
                r"Gs2 ($n_z=64$, $n_{period}=4$)",
                r"stella ($n_z=32$, $n_{period}=2$)",
                r"stella ($n_z=64$, $n_{period}=4$)",
                ],
                omega_diff=True,
                save_name="for_talk_cbc_Omegadiff_beta_scan.eps"

                                )
    return

def make_gbar_benchmark_for_thesis():
    """ """

    make_beta_plots_from_pickles_for_thesis_gbar(
        [
        "sims/" + mps.gs2_beta_scan_ky_05_np2_nt32_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + mps.gs2_beta_scan_ky_05_np4_nt64_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        ],
        [
        r"GS2 ($n_z=32$, $n_{period}=2$)",
        r"GS2 ($n_z=64$, $n_{period}=4$)",
        r"stella (implicit) ($n_z=32$, $n_{period}=2$)",
        r"stella (implicit) ($n_z=64$, $n_{period}=4$)",
        r"stella (implicit) ($n_z=64$, $n_{period}=2$)",
        r"stella (implicit) ($n_z=128$, $n_{period}=2$)",
        r"stella (implicit) ($n_z=256$, $n_{period}=2$)",
        r"stella (implicit) ($n_z=512$, $n_{period}=2$)",
        ],
        omega_diff=False,
        save_name="cbc_Omega_beta_scan_gbar.eps"

                        )

    make_beta_plots_from_pickles_for_thesis_gbar(
        [
        "sims/" + mps.gs2_beta_scan_ky_05_np2_nt32_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + mps.gs2_beta_scan_ky_05_np4_nt64_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        ],
        [
        r"GS2 ($n_z=32$, $n_{period}=2$)",
        r"GS2 ($n_z=64$, $n_{period}=4$)",
        r"stella (implicit) ($n_z=32$, $n_{period}=2$)",
        r"stella (implicit) ($n_z=64$, $n_{period}=4$)",
        r"stella (implicit) ($n_z=64$, $n_{period}=2$)",
        r"stella (implicit) ($n_z=128$, $n_{period}=2$)",
        r"stella (implicit) ($n_z=256$, $n_{period}=2$)",
        r"stella (implicit) ($n_z=512$, $n_{period}=2$)",
        ],
        omega_diff=True,
        save_name="cbc_Omega_beta_scan_gbar_diff.eps"

                        )
    return


def make_gbar_plot_for_talk():
    """ """

    make_beta_plots_from_pickles_for_thesis_gbar(
        [
        "sims/" + mps.gs2_beta_scan_ky_05_np2_nt32_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + mps.gs2_beta_scan_ky_05_np4_nt64_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        "sims/" + stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
        ],
        [
        r"GS2 ($n_z=32$, $n_{period}=2$)",
        r"GS2 ($n_z=64$, $n_{period}=4$)",
        r"stella ($n_z=32$, $n_{period}=2$)",
        r"stella ($n_z=64$, $n_{period}=4$)",
        r"stella ($n_z=64$, $n_{period}=2$)",
        r"stella ($n_z=128$, $n_{period}=2$)",
        r"stella ($n_z=256$, $n_{period}=2$)",
        r"stella ($n_z=512$, $n_{period}=2$)",
        ],
        omega_diff=False,
        save_name="for_talk_cbc_Omega_beta_scan_gbar.eps", for_talk=True

                        )


if __name__ == "__main__":
    print("Hello world")
    #make_omega_beta_plots_fapar1_fbpar1()
    make_gbar_benchmark_for_thesis()
    #make_gbar_plot_for_talk()
    #for_talk_make_omega_beta_plots_fapar1_fbpar1()
