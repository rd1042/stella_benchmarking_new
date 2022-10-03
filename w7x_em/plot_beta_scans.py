""" """

import numpy as np
import pickle
import matplotlib.pyplot as plt
from plot_2d_utils import uniquearrays2meshgrids
import math
import sys
import make_beta_scans as make_scans

# def beta_scan_gene_benchmark_for_thesis():
#     """ """
#
#     folder_longnames = [make_scans.folder_name_dkh_explicit_betalr,
#                 make_scans.folder_name_dkh_explicit_betalr_fbpar0,
#                 make_scans.folder_name_dkh_implicit_betalr,
#                 ]
#     labels=[r"explicit", r"explicit ($\delta B_\parallel=0$)", r"implicit"]
#
#     fig = plt.figure(figsize=(12,12))
#
#     left = 0.1
#     right = 0.96
#     top = 0.96
#     bottom = 0.1
#     vspace = 0.05
#     width = right - left
#     height = (top - bottom  - 2*vspace)/3
#     row2_bottom = bottom + height + vspace
#     row1_bottom = row2_bottom + height + vspace
#     ax1 = fig.add_axes((left, row1_bottom, width, height))
#     ax2 = fig.add_axes((left, row2_bottom, width, height))
#     ax3 = fig.add_axes((left, bottom, width, height))
#
#     label_fontsize = 40
#     legend_fontsize = 15
#     xticklabel_fontsize = 15
#     yticklabel_fontsize = 15
#     for counter, folder_longname in enumerate(folder_longnames):
#         pickle_longname = folder_longname + "/omega_beta_array.pickle"
#         file = open(pickle_longname, "rb")
#         [pickle_string, beta_array, gammaom_array, gammasafe_array,
#             freq_array] = pickle.load(file)
#         file.close()
#         ax1.plot(beta_array*100, -freq_array, marker="o")
#         ax2.plot(beta_array*100, gammaom_array, marker="o")
#         ax3.plot(beta_array*100, gammasafe_array, marker="o", label=labels[counter])
#
#     ax3.set_xlabel(r"$\beta$ (%)", fontsize=label_fontsize)
#     ax1.set_ylabel(r"$\tilde{\omega}$", fontsize=label_fontsize)
#     ax2.set_ylabel(r"$\tilde{\gamma}$", fontsize=label_fontsize)
#     ax3.set_ylabel(r"$\tilde{\gamma}_{2}$", fontsize=label_fontsize)
#     ax3.legend(loc="upper center", fontsize=legend_fontsize)
#     ax1.set_ylim(0, 0.5)
#     ax2.set_ylim(-0.18, 0.4)
#     ax3.set_ylim(0, 0.4)
#     for ax in [ax1, ax2, ax3]:
#         ax.set_xlim(-0.05, 20.05)
#         ax.set_xticks([0, 4, 8, 12, 16, 20])
#         ax.set_xticklabels([r"$0$", r"$4$", r"$8$", r"$12$", r"$16$", r"$20$"], fontsize=xticklabel_fontsize)
#
#     plt.savefig("w7x_gene_benchmark.eps")
#     plt.close()
#     return


def make_beta_scan_plot(folder_longnames):
    """For looking at the data, not fancy plotting """

    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)

    for folder_longname in folder_longnames:
        pickle_longname = folder_longname + "/omega_beta_array.pickle"
        file = open(pickle_longname, "rb")
        [pickle_string, beta_array, gammaom_array, gammasafe_array,
            freq_array] = pickle.load(file)
        file.close()
        ax1.plot(beta_array, freq_array, marker="o")
        ax2.plot(beta_array, gammaom_array, marker="o")
        ax3.plot(beta_array, gammasafe_array, marker="o")
    plt.show()
    return

def make_beta_scan_plot_for_thesis():
    """Want to scan beta for ky=1, 1.5, 2, 3 (or maybe fewer betas if time constrained)
    - For each ky, plot Omega(beta)
    - For some select ky, beta, plot:
        - phi(z)
        - apar(z)
        - bpar(z)
       together with B(z) """

    folder_longnames =  [
                make_scans.folder_name_w003_explicit,
                make_scans.folder_name_w003_implicit_ky1,
                make_scans.folder_name_w003_implicit_ky2,
                make_scans.folder_name_w003_implicit_ky15,
                make_scans.folder_name_w003_explicit_fbpar0,
                make_scans.folder_name_w003_implicit_ky1_fbpar0,
                make_scans.folder_name_w003_implicit_ky2_fbpar0,
                make_scans.folder_name_w003_implicit_ky15_fbpar0,
                        ]

    return

if __name__ == "__main__":
    print("Hello world")

    # make_beta_scan_plot([make_scans.folder_name_dkh_explicit_betalr,
    #                 make_scans.folder_name_dkh_explicit_betalr_fbpar0,
    #                 make_scans.folder_name_dkh_implicit_betalr,
    #                 ])
    # make_beta_scan_plot(make_scans.folder_name_dkh_explicit_betalr_fbpar0)
    # beta_scan_gene_benchmark_for_thesis()
    make_beta_scan_plot_for_thesis()
