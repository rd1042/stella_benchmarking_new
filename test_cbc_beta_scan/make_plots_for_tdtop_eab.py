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

def make_beta_plots_from_pickles_for_tdtop(pickle_longnames, labels, marker_size=1, omega_diff=False,
                                save_name="test.eps"):
    """Either plot Omega(beta), or, if omega_diff=True, plot (Omega-Omega_ref) (beta),
    where Omega_ref is taken from the first pickle. """

    marker_size = 12
    label_fontsize = 28
    legend_fontsize = 13
    marker_list = ["s", "o", "P", "X", "v", "^", "<", ">", "1", "2", "3"]
    lw_list = [4, 3, 3, 2]
    ls_list = ["-", "--", "-.", (0, (4,1,2,1))]

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

    top = 0.9
    left = 0.12
    right = 0.98
    bottom = 0.1
    vspace = 0.02
    hspace = 0.12
    height = (top - bottom - vspace)/2
    width = (right - left - hspace)/2
    row1_bottom = bottom + height + vspace
    col2_left = left + width + hspace

    print("col2_left, width = ", col2_left, width)
    fig = plt.figure(figsize=(15,8))
    ax1 = fig.add_axes((left, row1_bottom, width, height))
    ax2 = fig.add_axes((left, bottom, width, height))
    ax3 = fig.add_axes((col2_left, row1_bottom, width, height))
    ax4 = fig.add_axes((col2_left, bottom, width, height))

    ## Use the last folder as the "correct" one
    reference_pickle_longname = pickle_longnames[-1]
    myfile = open(reference_pickle_longname, "rb")
    [beta_ref, gamma_ref, omega_ref] = pickle.load(myfile)
    myfile.close()

    for folder_idx, pickle_longname in enumerate(pickle_longnames):
        myfile = open(pickle_longname, "rb")
        [beta_vals, gamma_vals, omega_vals] = pickle.load(myfile)
        myfile.close()
        ax1.plot(beta_vals, omega_vals, label=labels[folder_idx],
            lw=lw_list[folder_idx], ls=ls_list[folder_idx], marker=marker_list[folder_idx], mfc="none", markersize=marker_size)
        ax2.plot(beta_vals, gamma_vals, label=labels[folder_idx],
            lw=lw_list[folder_idx], ls=ls_list[folder_idx], marker=marker_list[folder_idx], mfc="none", markersize=marker_size)
        if folder_idx!=(len(pickle_longnames)-1):
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
            ax3.plot(beta_to_compare, omega_diff, label=labels[folder_idx],
                     lw=lw_list[folder_idx], ls=ls_list[folder_idx], marker=marker_list[folder_idx], mfc="none", markersize=marker_size)
            ax4.plot(beta_to_compare, gamma_diff, label=labels[folder_idx],
                     lw=lw_list[folder_idx], ls=ls_list[folder_idx], marker=marker_list[folder_idx], mfc="none", markersize=marker_size)

    # for ax in [ax1, ax2]:
    #     ax.grid(True)


    ax2.legend(loc="best", fontsize=legend_fontsize)
    min_beta = -0.001
    max_beta = 0.041
    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_xlim((min_beta, max_beta))
        ax.tick_params("y", which="major", length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.set_xticks([0, 0.01, 0.02, 0.03, 0.04])
        ax.set_xticklabels([r"$0$", r"$0.01$", r"$0.02$", r"$0.03$", r"$0.04$"], fontsize=x_ticklabelfontsize)
        ax.tick_params("x", which="major", length=my_xticklength, width=my_xtickwidth, direction="in", labelbottom=True)
    for ax in [ax1, ax3]:
        ax.tick_params("x", which="major", length=my_xticklength, width=my_xtickwidth, direction="in", labelbottom=False)

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
    ax1.set_ylabel(r"frequency" + "\n" + r"$(v_{th,i}/a)$", fontsize=label_fontsize)
    ax2.set_ylabel(r"growth rate" + "\n" + r"$(v_{th,i}/a)$", fontsize=label_fontsize)
    # ax1.annotate(r"$\{$", (0.002, 0.4), rotation=90, fontsize=bracket_fontsize)
    itg_x = 0.005
    kbm_x = 0.027
    ax1.annotate(r"ITG", xy=(itg_x, 0.4), xytext=(itg_x, 0.43), fontsize=itg_fontsize,
                ha='center', va='bottom',
                    arrowprops=dict(arrowstyle='-[, widthB=1.8, lengthB=0.5', lw=2.0) )
    ax1.annotate(r"KBM", xy=(kbm_x, 0.38), xytext=(kbm_x, 0.35), fontsize=itg_fontsize,
                ha='center', va='top',
                    arrowprops=dict(arrowstyle='-[, widthB=4.0, lengthB=0.5', lw=2.0) )


    ax3.plot((min_beta, max_beta), (0,0), c="gray", zorder=-10)
    ax4.plot((min_beta, max_beta), (0,0), c="gray", zorder=-10)
    ax3.set_title(r"comparison with " + "\n" + r"GS2 ($n_\theta=64$, $n_{p}=4$)", fontsize=20)
    ax3.set_ylabel(r"$\%$ diff." + "\n" + r"(freq.)", fontsize=label_fontsize)
    ax4.set_ylabel(r"$\%$ diff." + "\n" + r"(growth rate)", fontsize=label_fontsize)
    ax3.set_yticks([0, 2, 4, 6, 8])
    ax3.set_yticklabels([r"$0$", r"$2$", r"$4$", r"$6$", r"$8$"], fontsize=x_ticklabelfontsize)
    ax4.set_yticks([-20, -10, 0, 10])
    ax4.set_yticklabels([r"$-20$", r"$-10$", r"$0$", r"$10$"], fontsize=x_ticklabelfontsize)
    ax4.set_xlabel(r"$\beta$", fontsize=label_fontsize)
    ax2.set_xlabel(r"$\beta$", fontsize=label_fontsize)
    plt.savefig(save_name)
    plt.close()

    return

def make_omega_beta_plots_fapar1_fbpar1_1():
    """ """

    make_beta_plots_from_pickles_for_tdtop(
                [
                "sims/" + mps.stella_src_h_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_zupw0_fapar1_fbpar1_str_mirror_implicit_folder + "/beta_gamma_omega.pickle",
                "sims/" + mps.stella_src_h_beta_scan_ky_05_np4_nt64_nvpa18_nmu12_zupw0_fapar1_fbpar1_str_mirror_implicit_folder + "/beta_gamma_omega.pickle",
                "sims/" + mps.gs2_beta_scan_ky_05_np2_nt32_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                "sims/" + mps.gs2_beta_scan_ky_05_np4_nt64_ng8_ne18_fapar1_fbpar1_folder + "/beta_gamma_omega.pickle",
                ],
                [
                r"stella ($n_\theta=32$, $n_{period}=2$)",
                r"stella ($n_\theta=64$, $n_{period}=4$)",
                r"GS2 ($n_\theta=32$, $n_{period}=2$)",
                r"GS2 ($n_\theta=64$, $n_{period}=4$)",
                ],
                omega_diff=False,
                save_name="cbc_Omega_beta_scan_tdotp.eps"

                                )
    return

if __name__ == "__main__":
    print("Hello world")
    make_omega_beta_plots_fapar1_fbpar1_1()
