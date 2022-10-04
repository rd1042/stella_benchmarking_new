""" """

import numpy as np
import pickle
import matplotlib.pyplot as plt
from plot_2d_utils import uniquearrays2meshgrids
import math
import sys
import make_beta_scans as make_scans
from save_pickles_from_stella_scans import get_omega_data, get_phiz_data_stella

sys.path.append("../postprocessing_tools")
import helper_linear_sims as help_lin
from extract_sim_data import get_omega_data, get_phiz_data, get_aparz_data, get_bparz_data
from helper_ncdf_new import view_ncdf_variables_with_xarray, extract_data_from_ncdf_with_xarray

default_cmap = plt.get_cmap("tab10")

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

def kjm3_make_beta_scan_plot_for_thesis():
    """Want to scan beta for ky=1, 1.5, 2, 3 (or maybe fewer betas if time constrained)
    - For each ky, plot Omega(beta)
    - For some select ky, beta, plot:
        - phi(z)
        - apar(z)
        - bpar(z)
       together with B(z) """

    folder_longnames =  [
                make_scans.folder_name_kjm3_explicit,
                make_scans.folder_name_kjm3_implicit_ky1_upw0,
                make_scans.folder_name_kjm3_implicit_ky15_upw0,
                make_scans.folder_name_kjm3_implicit_ky2_upw0,
                make_scans.folder_name_kjm3_explicit_fbpar0,
                make_scans.folder_name_kjm3_implicit_ky1_fbpar0_upw0,
                        ]
    labels = [ r"$\tilde{k}_y=1$, RK2",
               r"$\tilde{k}_y=1$, implicit",
               r"$\tilde{k}_y=1.5$, implicit",
               r"$\tilde{k}_y=2$, implicit",
               r"$\tilde{k}_y=1$, $B_\parallel=0$, RK2",
               r"$\tilde{k}_y=1$, $B_\parallel=0$, implicit",

                ]
    special_cols = [default_cmap(6), default_cmap(7),
                default_cmap(8), default_cmap(9),
                ]
    special_beta_idxs = [4, 6, 8, 10]
    fig = plt.figure(figsize=(12,12))
    left = 0.1
    right = 0.98
    bottom = 0.1
    top = 0.985
    vspace = 0.015
    width = right - left
    height = (top - bottom - 2*vspace)/3

    row2_bottom = bottom + height + vspace
    row1_bottom = row2_bottom + height + vspace

    x_ticklabelfontsize = 20
    x_labelfontsize = 40
    my_xticklength = 6
    my_xtickwidth = 3
    my_xticklength_minor = 3
    my_xtickwidth_minor = 1
    y_labelfontsize = x_labelfontsize
    y_ticklabelfontsize = x_ticklabelfontsize
    # cb_ticklabelsize = 14
    # cbax_title_fontsize = 20
    ax_title_fontsize = 20
    marker_list = ["s", "o", "P", "X", "v", "^", "<", ">", "1", "2", "3"]
    marker_size=15
    my_linewidth = 3
    legend_fontsize = 13.5

    ax1 = fig.add_axes((left, row1_bottom, width, height))
    ax2 = fig.add_axes((left, row2_bottom, width, height))
    ax3 = fig.add_axes((left, bottom, width, height))

    for folder_idx, folder_longname in enumerate(folder_longnames):
        pickle_longname = folder_longname + "/omega_beta_array.pickle"
        file = open(pickle_longname, "rb")
        [pickle_string, beta_array, gammaom_array, gammasafe_array,
            freq_array] = pickle.load(file)
        file.close()
        ax1.plot(beta_array, -freq_array, marker=marker_list[folder_idx], lw=my_linewidth,
                    markersize=marker_size, mfc="None", label=labels[folder_idx], c=default_cmap(folder_idx))
        ax2.plot(beta_array, gammaom_array, marker=marker_list[folder_idx], lw=my_linewidth,
                    markersize=marker_size, mfc="None", label=labels[folder_idx], c=default_cmap(folder_idx))
        ax3.plot(beta_array, gammasafe_array, marker=marker_list[folder_idx], lw=my_linewidth,
                    markersize=marker_size, mfc="None", label=labels[folder_idx], c=default_cmap(folder_idx))
        if folder_idx == 1:
            for counter, beta_idx in enumerate(special_beta_idxs):
                ax1.plot(beta_array[beta_idx], -freq_array[beta_idx], marker=marker_list[folder_idx], lw=my_linewidth,
                    markersize=marker_size, mfc=special_cols[counter], c=default_cmap(folder_idx))
                ax2.plot(beta_array[beta_idx], gammaom_array[beta_idx], marker=marker_list[folder_idx], lw=my_linewidth,
                    markersize=marker_size, mfc=special_cols[counter], c=default_cmap(folder_idx))
                ax3.plot(beta_array[beta_idx], gammasafe_array[beta_idx], marker=marker_list[folder_idx], lw=my_linewidth,
                    markersize=marker_size, mfc=special_cols[counter], c=default_cmap(folder_idx))

    for ax in [ax1,ax2,ax3]:
        ax.set_xlim((-0.001,0.101))
    ax1.set_ylim(0, 2.8)
    ax2.set_ylim(-0.3, 3.9)
    ax3.set_ylim(-0.3, 3.9)

    for ax in [ax1, ax2, ax3]:
        #ax.scatter(fprim_list, tprim_list, c="r", s=3.)
        #print("in make_plots. sorted(set(fprim_list)) = ", sorted(set(fprim_list)))
        ax.tick_params("x", which="major", top=True, bottom=True, length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.tick_params("y", which="major", left=True, right=True, length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.tick_params("x", which="minor", top=True, bottom=True, length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in")
        ax.tick_params("y", which="minor", left=True, right=True, length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in")
        ax.set_xticks([0, 0.02, 0.04, 0.06, 0.08, 0.1])
        ax.set_xticks([0.01, 0.03, 0.05, 0.07, 0.09], minor=True)
        ax.set_xticklabels([], minor=True)
    ax1.set_yticks([0, 1, 2])
    ax1.set_yticklabels([r"$0$", r"$1$", r"$2$"], fontsize=y_ticklabelfontsize)
    ax1.set_yticks([0.5, 1.5, 2.5], minor=True)
    ax1.set_yticklabels([], minor=True)

    for ax in [ax2, ax3]:
        ax.set_yticks([0, 1, 2, 3])
        ax.set_yticklabels([r"$0$", r"$1$", r"$2$", r"$3$"], fontsize=y_ticklabelfontsize)
        ax.set_yticks([0.5, 1.5, 2.5, 3.5], minor=True)
        ax.set_yticklabels([], minor=True)

    ax3.legend(loc="upper left", fontsize=legend_fontsize)
    ax3.set_xticklabels([r"$0$", r"$2$", r"$4$", r"$6$", r"$8$", r"$10$"], fontsize=x_ticklabelfontsize)
    ax1.set_xticklabels([], fontsize=x_ticklabelfontsize)
    ax2.set_xticklabels([], fontsize=x_ticklabelfontsize)
    ax3.set_xlabel(r"$\beta (\%)$", fontsize=x_labelfontsize)
    ax1.set_ylabel(r"$\tilde{\omega}_{stella}$", fontsize=y_labelfontsize)
    ax2.set_ylabel(r"$\tilde{\gamma}_{stella}$", fontsize=y_labelfontsize)
    ax3.set_ylabel(r"$\tilde{\gamma}_{2}$", fontsize=y_labelfontsize)

    plt.savefig("images/kjm3_omega_beta_scan.eps")
    plt.close()

    return

def w003_make_beta_scan_plot_for_thesis():
    """Want to scan beta for ky=1, 1.5, 2, 3 (or maybe fewer betas if time constrained)
    - For each ky, plot Omega(beta)
    - For some select ky, beta, plot:
        - phi(z)
        - apar(z)
        - bpar(z)
       together with B(z) """

    folder_longnames =  [
                make_scans.folder_name_w003_explicit,
                make_scans.folder_name_w003_implicit_ky1_upw0,
                make_scans.folder_name_w003_implicit_ky15_upw0,
                make_scans.folder_name_w003_implicit_ky2_upw0,
                make_scans.folder_name_w003_explicit_fbpar0,
                make_scans.folder_name_w003_implicit_ky1_fbpar0_upw0,
                        ]
    labels = [ r"$\tilde{k}_y=1$, RK2",
               r"$\tilde{k}_y=1$, implicit",
               r"$\tilde{k}_y=1.5$, implicit",
               r"$\tilde{k}_y=2$, implicit",
               r"$\tilde{k}_y=1$, $B_\parallel=0$, RK2",
               r"$\tilde{k}_y=1$, $B_\parallel=0$, implicit",

                ]
    special_cols = [default_cmap(6), default_cmap(7),
                default_cmap(8), default_cmap(9),
                ]
    special_beta_idxs = [4, 6, 8, 10]
    fig = plt.figure(figsize=(12,12))
    left = 0.1
    right = 0.98
    bottom = 0.1
    top = 0.985
    vspace = 0.015
    width = right - left
    height = (top - bottom - 2*vspace)/3

    row2_bottom = bottom + height + vspace
    row1_bottom = row2_bottom + height + vspace

    x_ticklabelfontsize = 20
    x_labelfontsize = 40
    my_xticklength = 6
    my_xtickwidth = 3
    my_xticklength_minor = 3
    my_xtickwidth_minor = 1
    y_labelfontsize = x_labelfontsize
    y_ticklabelfontsize = x_ticklabelfontsize
    # cb_ticklabelsize = 14
    # cbax_title_fontsize = 20
    ax_title_fontsize = 20
    marker_list = ["s", "o", "P", "X", "v", "^", "<", ">", "1", "2", "3"]
    marker_size=15
    my_linewidth = 3
    legend_fontsize = 13.5

    ax1 = fig.add_axes((left, row1_bottom, width, height))
    ax2 = fig.add_axes((left, row2_bottom, width, height))
    ax3 = fig.add_axes((left, bottom, width, height))

    for folder_idx, folder_longname in enumerate(folder_longnames):
        pickle_longname = folder_longname + "/omega_beta_array.pickle"
        file = open(pickle_longname, "rb")
        [pickle_string, beta_array, gammaom_array, gammasafe_array,
            freq_array] = pickle.load(file)
        file.close()
        ax1.plot(beta_array, -freq_array, marker=marker_list[folder_idx], lw=my_linewidth,
                    markersize=marker_size, mfc="None", label=labels[folder_idx], c=default_cmap(folder_idx))
        ax2.plot(beta_array, gammaom_array, marker=marker_list[folder_idx], lw=my_linewidth,
                    markersize=marker_size, mfc="None", label=labels[folder_idx], c=default_cmap(folder_idx))
        ax3.plot(beta_array, gammasafe_array, marker=marker_list[folder_idx], lw=my_linewidth,
                    markersize=marker_size, mfc="None", label=labels[folder_idx], c=default_cmap(folder_idx))
        if folder_idx == 1:
            for counter, beta_idx in enumerate(special_beta_idxs):
                ax1.plot(beta_array[beta_idx], -freq_array[beta_idx], marker=marker_list[folder_idx], lw=my_linewidth,
                    markersize=marker_size, mfc=special_cols[counter], c=default_cmap(folder_idx))
                ax2.plot(beta_array[beta_idx], gammaom_array[beta_idx], marker=marker_list[folder_idx], lw=my_linewidth,
                    markersize=marker_size, mfc=special_cols[counter], c=default_cmap(folder_idx))
                ax3.plot(beta_array[beta_idx], gammasafe_array[beta_idx], marker=marker_list[folder_idx], lw=my_linewidth,
                    markersize=marker_size, mfc=special_cols[counter], c=default_cmap(folder_idx))

    for ax in [ax1,ax2,ax3]:
        ax.set_xlim((-0.001,0.101))
    ax1.set_ylim(0, 2.8)
    ax2.set_ylim(-0.3, 3.9)
    ax3.set_ylim(-0.3, 3.9)

    for ax in [ax1, ax2, ax3]:
        #ax.scatter(fprim_list, tprim_list, c="r", s=3.)
        #print("in make_plots. sorted(set(fprim_list)) = ", sorted(set(fprim_list)))
        ax.tick_params("x", which="major", top=True, bottom=True, length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.tick_params("y", which="major", left=True, right=True, length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.tick_params("x", which="minor", top=True, bottom=True, length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in")
        ax.tick_params("y", which="minor", left=True, right=True, length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in")
        ax.set_xticks([0, 0.02, 0.04, 0.06, 0.08, 0.1])
        ax.set_xticks([0.01, 0.03, 0.05, 0.07, 0.09], minor=True)
        ax.set_xticklabels([], minor=True)
    ax1.set_yticks([0, 1, 2])
    ax1.set_yticklabels([r"$0$", r"$1$", r"$2$"], fontsize=y_ticklabelfontsize)
    ax1.set_yticks([0.5, 1.5, 2.5], minor=True)
    ax1.set_yticklabels([], minor=True)

    for ax in [ax2, ax3]:
        ax.set_yticks([0, 1, 2, 3])
        ax.set_yticklabels([r"$0$", r"$1$", r"$2$", r"$3$"], fontsize=y_ticklabelfontsize)
        ax.set_yticks([0.5, 1.5, 2.5, 3.5], minor=True)
        ax.set_yticklabels([], minor=True)

    ax3.legend(loc="upper left", fontsize=legend_fontsize)
    ax3.set_xticklabels([r"$0$", r"$2$", r"$4$", r"$6$", r"$8$", r"$10$"], fontsize=x_ticklabelfontsize)
    ax1.set_xticklabels([], fontsize=x_ticklabelfontsize)
    ax2.set_xticklabels([], fontsize=x_ticklabelfontsize)
    ax3.set_xlabel(r"$\beta (\%)$", fontsize=x_labelfontsize)
    ax1.set_ylabel(r"$\tilde{\omega}_{stella}$", fontsize=y_labelfontsize)
    ax2.set_ylabel(r"$\tilde{\gamma}_{stella}$", fontsize=y_labelfontsize)
    ax3.set_ylabel(r"$\tilde{\gamma}_{2}$", fontsize=y_labelfontsize)

    plt.savefig("images/w003_omega_beta_scan.eps")
    plt.close()

    return

def get_normalised_fields_for_z_zeta(sim_longname):
    """ """
    z, real_phi, imag_phi = get_phiz_data(sim_longname, "stella")
    abs_phi = help_lin.get_abs(real_phi, imag_phi)
    z, real_apar, imag_apar = get_aparz_data(sim_longname, "stella")
    abs_apar = help_lin.get_abs(real_apar, imag_apar)
    z, real_bpar, imag_bpar = get_bparz_data(sim_longname, "stella")
    abs_bpar = help_lin.get_abs(real_bpar, imag_bpar)
    max_abs_phi = np.max(abs_phi)
    abs_phi = abs_phi/max_abs_phi
    abs_apar = abs_apar/max_abs_phi
    abs_bpar = abs_bpar/max_abs_phi

    geo_longname = sim_longname + ".geometry"
    data = np.loadtxt(geo_longname, skiprows=4, usecols=(1,2,3)) # data = [zed, zeta, bmaf]
    zed = data[:,0]
    zeta = data[:,1]
    abs_phi = np.sqrt(real_phi*real_phi + imag_phi*imag_phi)
    max_abs_phi = np.max(abs_phi)
    real_phi = real_phi/max_abs_phi
    imag_phi = imag_phi/max_abs_phi
    abs_phi = abs_phi/max_abs_phi

    return z, zeta, abs_phi, abs_apar, abs_bpar

def kjm3_make_mode_plots_for_beta_scan_for_thesis():
    """Want to scan beta for ky=1, 1.5, 2, 3 (or maybe fewer betas if time constrained)
    - For each ky, plot Omega(beta)
    - For some select ky, beta, plot:
        - phi(z)
        - apar(z)
        - bpar(z)
       together with B(z) """

    sim_longnames =  [
                # make_scans.folder_name_kjm3_implicit_ky1_upw0 + "/beta_0.0400/input",
                # make_scans.folder_name_kjm3_implicit_ky1_upw0 + "/beta_0.0600/input",
                # make_scans.folder_name_kjm3_implicit_ky1_upw0 + "/beta_0.0800/input",
                # make_scans.folder_name_kjm3_implicit_ky1_upw0 + "/beta_0.1000/input",
                make_scans.folder_name_kjm3_implicit_ky1_upw0 + "/beta_0.0400/input_for_fields",
                make_scans.folder_name_kjm3_implicit_ky1_upw0 + "/beta_0.0600/input_for_fields",
                make_scans.folder_name_kjm3_implicit_ky1_upw0 + "/beta_0.0800/input_for_fields",
                make_scans.folder_name_kjm3_implicit_ky1_upw0 + "/beta_0.1000/input_for_fields",
                #make_scans.folder_name_kjm3_implicit_ky1_upw0 + "/beta_0.1000/input_different_res",
                        ]
    labels = [ r"$\beta=4\%$",
               r"$\beta=6\%$",
               r"$\beta=8\%$",
               r"$\beta=10\%$",
               #r"$\beta=10\%$ (new resolution)",
                ]
    cols = [default_cmap(6), default_cmap(7),
            default_cmap(8), default_cmap(9),
            ]

    lw_list = [6, 4, 4, 2, 2]
    ls_list = ["-", "--", "-.", (0, (4,1,2,1)), "-"]

    fig = plt.figure(figsize=(12,12))
    left = 0.14
    right = 0.98
    bottom = 0.08
    top = 0.9
    vspace = 0.015
    width = right - left
    height = (top - bottom - 2*vspace)/3

    row2_bottom = bottom + height + vspace
    row1_bottom = row2_bottom + height + vspace

    x_ticklabelfontsize = 20
    x_labelfontsize = 40
    my_xticklength = 6
    my_xtickwidth = 3
    my_xticklength_minor = 3
    my_xtickwidth_minor = 1
    y_labelfontsize = x_labelfontsize
    y_ticklabelfontsize = x_ticklabelfontsize
    inset_ticklabelfonstsize = 16
    # cb_ticklabelsize = 14
    # cbax_title_fontsize = 20
    ax_title_fontsize = 20
    marker_list = ["s", "o", "P", "X", "v", "^", "<", ">", "1", "2", "3"]
    marker_size=15
    my_linewidth = 3
    legend_fontsize = 20

    inset_width = 0.26
    inset_height = 0.16
    inset_left = 0.04
    inset_bottom = 0.08

    ax1 = fig.add_axes((left, row1_bottom, width, height))
    ax2 = fig.add_axes((left, row2_bottom, width, height))
    ax3 = fig.add_axes((left, bottom, width, height))
    ax1_inset = fig.add_axes((left + inset_left, row1_bottom + inset_bottom, inset_width, inset_height))
    ax1zeta = ax1.twiny() # fig.add_axes((left, row1_bottom, width, height))
    ax2zeta = ax2.twiny() # fig.add_axes((left, row1_bottom, width, height))
    ax3zeta = ax3.twiny() # fig.add_axes((left, row1_bottom, width, height))

    for sim_idx, sim_longname in enumerate(sim_longnames):
        z, zeta, abs_phi, abs_apar, abs_bpar = get_normalised_fields_for_z_zeta(sim_longname)
        if sim_idx == 4:
            z = z * 0.2
        if sim_idx == 0:
            ax1_inset.plot(z/np.pi, abs_phi, lw=lw_list[sim_idx], ls=ls_list[sim_idx],
                        label=labels[sim_idx], c=cols[sim_idx])
        # print("abs_phi = ", abs_phi)
        # print("min(zeta/pi) = ", np.min(zeta/np.pi))
        # print("max(zeta/pi) = ", np.max(zeta/np.pi))
        ax1.plot(z/np.pi, abs_phi, lw=lw_list[sim_idx], ls=ls_list[sim_idx],
                    label=labels[sim_idx], c=cols[sim_idx])
        ax2.plot(z/np.pi, abs_apar, lw=lw_list[sim_idx], ls=ls_list[sim_idx],
                    label=labels[sim_idx], c=cols[sim_idx])
        ax3.plot(z/np.pi, abs_bpar, lw=lw_list[sim_idx], ls=ls_list[sim_idx],
                    label=labels[sim_idx], c=cols[sim_idx])

    # for ax in [ax1,ax2,ax3]:
    #     ax.set_xlim((0,0.1))
    # ax1.set_ylim(0, 2.8)
    # ax2.set_ylim(-0.3, 3.9)
    # ax3.set_ylim(-0.3, 3.9)

    for ax in [ax1, ax2, ax3]:
        #ax.scatter(fprim_list, tprim_list, c="r", s=3.)
        #print("in make_plots. sorted(set(fprim_list)) = ", sorted(set(fprim_list)))
        ax.set_xlim(-0.15, 0.15)
        ax.tick_params("x", which="major", top=False, bottom=True, length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.tick_params("y", which="major", left=True, right=True, length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.tick_params("x", which="minor", top=False, bottom=True, length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in")
        ax.tick_params("y", which="minor", left=True, right=True, length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in")
        ax.set_xticks([-0.1, 0, 0.1])
        ax.set_xticks([-0.15, -0.05, 0.05, 0.15], minor=True)
        ax.set_xticklabels([], minor=True)
        ax.grid(True)

    ax3.set_xlabel(r"$z$", fontsize=x_labelfontsize)
    ax1zeta.set_xlabel(r"$\zeta$", fontsize=x_labelfontsize, labelpad=20)
    #ax1zeta.set_xlabel(r"$\zeta$", fontsize=x_labelfontsize)
    ax1.set_ylabel(r"$\vert \tilde{\phi}_k \vert$", fontsize=y_labelfontsize)
    ax2.set_ylabel(r"$\vert \tilde{A}_{\parallel k} \vert$", fontsize=y_labelfontsize)
    ax3.set_ylabel(r"$\vert \tilde{B}_{\parallel k} \vert$", fontsize=y_labelfontsize)
    ax1.set_yticks([0, 0.5, 1,])
    ax1.set_yticklabels([r"$0$", r"$0.5$", r"$1$"], fontsize=y_ticklabelfontsize)
    ax1.set_yticks([0.25, 0.75], minor=True)

    ax2.set_yticks([0, 0.05])
    ax2.set_yticklabels([r"$0$", r"$0.05$"], fontsize=y_ticklabelfontsize)
    ax2.set_yticks([0.25], minor=True)
    ax2.set_yticklabels([], minor=True)
    ax3.set_yticks([0, 0.05])
    ax3.set_yticklabels([r"$0$", r"$0.05$"], fontsize=y_ticklabelfontsize)
    ax3.set_yticks([0.25], minor=True)
    ax3.set_yticklabels([], minor=True)
    ax3.set_xticklabels([r"$-0.1$", "$0$", r"$0.1$"], fontsize=x_ticklabelfontsize)
    ax1.set_xticklabels([], fontsize=x_ticklabelfontsize)
    ax2.set_xticklabels([], fontsize=x_ticklabelfontsize)

    for ax in [ax1zeta, ax2zeta, ax3zeta]:
        ax.tick_params("x", which="major", top=True, bottom=False, length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.tick_params("x", which="minor", top=True, bottom=False, length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in")
        ax.set_xlim(np.min(zeta)*0.15/np.pi, np.max(zeta)*0.15/np.pi)
        ax.set_xticks([-0.4, 0, 0.4])
        ax.set_xticks([-0.6, -0.2, 0.2, 0.6], minor=True)
    ax1zeta.set_xticklabels([r"$-2\pi/5$", r"$0$", r"$2\pi/5$"], fontsize=x_ticklabelfontsize)
    ax2zeta.set_xticklabels([], fontsize=x_ticklabelfontsize)
    ax3zeta.set_xticklabels([], fontsize=x_ticklabelfontsize)
    # for ax in [ax2, ax3]:
    #     ax.set_yticks([0, 1, 2, 3])
    #     ax.set_yticklabels([r"$0$", r"$1$", r"$2$", r"$3$"], fontsize=y_ticklabelfontsize)
    #     ax.set_yticks([0.5, 1.5, 2.5, 3.5], minor=True)
    #     ax.set_yticklabels([], minor=True)

    ax3.legend(loc="upper left", fontsize=legend_fontsize)
    # ax3.set_xticklabels([r"$0$", r"$2$", r"$4$", r"$6$", r"$8$", r"$10$"], fontsize=x_ticklabelfontsize)
    # ax1.set_xticklabels([], fontsize=x_ticklabelfontsize)
    # ax2.set_xticklabels([], fontsize=x_ticklabelfontsize)

    ax1_inset.set_xlim((-1,1))
    ax1_inset.set_ylim((-0.05, 1.05))
    ax1_inset.set_yticks([0,1])
    ax1_inset.set_yticklabels([r"$0$",r"$1$"], fontsize=inset_ticklabelfonstsize)
    ax1_inset.set_xticks([-1,0,1])
    ax1_inset.set_xticklabels([r"$-1$", r"$0$",r"$1$"], fontsize=inset_ticklabelfonstsize)
    plt.savefig("images/kjm3_fields_beta_scan.eps")
    plt.close()

    return

def w003_make_mode_plots_for_beta_scan_for_thesis():
    """Want to scan beta for ky=1, 1.5, 2, 3 (or maybe fewer betas if time constrained)
    - For each ky, plot Omega(beta)
    - For some select ky, beta, plot:
        - phi(z)
        - apar(z)
        - bpar(z)
       together with B(z) """

    sim_longnames =  [
                # make_scans.folder_name_kjm3_implicit_ky1_upw0 + "/beta_0.0400/input",
                # make_scans.folder_name_kjm3_implicit_ky1_upw0 + "/beta_0.0600/input",
                # make_scans.folder_name_kjm3_implicit_ky1_upw0 + "/beta_0.0800/input",
                # make_scans.folder_name_kjm3_implicit_ky1_upw0 + "/beta_0.1000/input",
                make_scans.folder_name_w003_implicit_ky1_upw0 + "/beta_0.0400/input_for_fields",
                make_scans.folder_name_w003_implicit_ky1_upw0 + "/beta_0.0600/input_for_fields",
                make_scans.folder_name_w003_implicit_ky1_upw0 + "/beta_0.0800/input_for_fields",
                make_scans.folder_name_w003_implicit_ky1_upw0 + "/beta_0.1000/input_for_fields",
                #make_scans.folder_name_kjm3_implicit_ky1_upw0 + "/beta_0.1000/input_different_res",
                        ]
    labels = [ r"$\beta=4\%$",
               r"$\beta=6\%$",
               r"$\beta=8\%$",
               r"$\beta=10\%$",
               #r"$\beta=10\%$ (new resolution)",
                ]
    cols = [default_cmap(6), default_cmap(7),
            default_cmap(8), default_cmap(9),
            ]

    lw_list = [6, 4, 4, 2, 2]
    ls_list = ["-", "--", "-.", (0, (4,1,2,1)), "-"]

    fig = plt.figure(figsize=(12,12))
    left = 0.14
    right = 0.98
    bottom = 0.08
    top = 0.9
    vspace = 0.015
    width = right - left
    height = (top - bottom - 2*vspace)/3

    row2_bottom = bottom + height + vspace
    row1_bottom = row2_bottom + height + vspace

    x_ticklabelfontsize = 20
    x_labelfontsize = 40
    my_xticklength = 6
    my_xtickwidth = 3
    my_xticklength_minor = 3
    my_xtickwidth_minor = 1
    y_labelfontsize = x_labelfontsize
    y_ticklabelfontsize = x_ticklabelfontsize
    inset_ticklabelfonstsize = 16
    # cb_ticklabelsize = 14
    # cbax_title_fontsize = 20
    ax_title_fontsize = 20
    marker_list = ["s", "o", "P", "X", "v", "^", "<", ">", "1", "2", "3"]
    marker_size=15
    my_linewidth = 3
    legend_fontsize = 20

    inset_width = 0.26
    inset_height = 0.16
    inset_left = 0.04
    inset_bottom = 0.08

    ax1 = fig.add_axes((left, row1_bottom, width, height))
    ax2 = fig.add_axes((left, row2_bottom, width, height))
    ax3 = fig.add_axes((left, bottom, width, height))
    ax1_inset = fig.add_axes((left + inset_left, row1_bottom + inset_bottom, inset_width, inset_height))
    ax1zeta = ax1.twiny() # fig.add_axes((left, row1_bottom, width, height))
    ax2zeta = ax2.twiny() # fig.add_axes((left, row1_bottom, width, height))
    ax3zeta = ax3.twiny() # fig.add_axes((left, row1_bottom, width, height))

    for sim_idx, sim_longname in enumerate(sim_longnames):
        z, zeta, abs_phi, abs_apar, abs_bpar = get_normalised_fields_for_z_zeta(sim_longname)
        if sim_idx == 4:
            z = z * 0.2
        if sim_idx == 0:
            ax1_inset.plot(z/np.pi, abs_phi, lw=lw_list[sim_idx], ls=ls_list[sim_idx],
                        label=labels[sim_idx], c=cols[sim_idx])
        # print("abs_phi = ", abs_phi)
        # print("min(zeta/pi) = ", np.min(zeta/np.pi))
        # print("max(zeta/pi) = ", np.max(zeta/np.pi))
        ax1.plot(z/np.pi, abs_phi, lw=lw_list[sim_idx], ls=ls_list[sim_idx],
                    label=labels[sim_idx], c=cols[sim_idx])
        ax2.plot(z/np.pi, abs_apar, lw=lw_list[sim_idx], ls=ls_list[sim_idx],
                    label=labels[sim_idx], c=cols[sim_idx])
        ax3.plot(z/np.pi, abs_bpar, lw=lw_list[sim_idx], ls=ls_list[sim_idx],
                    label=labels[sim_idx], c=cols[sim_idx])

    # for ax in [ax1,ax2,ax3]:
    #     ax.set_xlim((0,0.1))
    # ax1.set_ylim(0, 2.8)
    # ax2.set_ylim(-0.3, 3.9)
    # ax3.set_ylim(-0.3, 3.9)

    for ax in [ax1, ax2, ax3]:
        #ax.scatter(fprim_list, tprim_list, c="r", s=3.)
        #print("in make_plots. sorted(set(fprim_list)) = ", sorted(set(fprim_list)))
        ax.set_xlim(-0.15, 0.15)
        ax.tick_params("x", which="major", top=False, bottom=True, length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.tick_params("y", which="major", left=True, right=True, length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.tick_params("x", which="minor", top=False, bottom=True, length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in")
        ax.tick_params("y", which="minor", left=True, right=True, length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in")
        ax.set_xticks([-0.1, 0, 0.1])
        ax.set_xticks([-0.15, -0.05, 0.05, 0.15], minor=True)
        ax.set_xticklabels([], minor=True)
        ax.grid(True)

    ax3.set_xlabel(r"$z$", fontsize=x_labelfontsize)
    ax1zeta.set_xlabel(r"$\zeta$", fontsize=x_labelfontsize, labelpad=20)
    #ax1zeta.set_xlabel(r"$\zeta$", fontsize=x_labelfontsize)
    ax1.set_ylabel(r"$\vert \tilde{\phi}_k \vert$", fontsize=y_labelfontsize)
    ax2.set_ylabel(r"$\vert \tilde{A}_{\parallel k} \vert$", fontsize=y_labelfontsize)
    ax3.set_ylabel(r"$\vert \tilde{B}_{\parallel k} \vert$", fontsize=y_labelfontsize)
    ax1.set_yticks([0, 0.5, 1,])
    ax1.set_yticklabels([r"$0$", r"$0.5$", r"$1$"], fontsize=y_ticklabelfontsize)
    ax1.set_yticks([0.25, 0.75], minor=True)

    ax2.set_yticks([0, 0.05])
    ax2.set_yticklabels([r"$0$", r"$0.05$"], fontsize=y_ticklabelfontsize)
    ax2.set_yticks([0.25], minor=True)
    ax2.set_yticklabels([], minor=True)
    ax3.set_yticks([0, 0.05])
    ax3.set_yticklabels([r"$0$", r"$0.05$"], fontsize=y_ticklabelfontsize)
    ax3.set_yticks([0.25], minor=True)
    ax3.set_yticklabels([], minor=True)
    ax3.set_xticklabels([r"$-0.1$", "$0$", r"$0.1$"], fontsize=x_ticklabelfontsize)
    ax1.set_xticklabels([], fontsize=x_ticklabelfontsize)
    ax2.set_xticklabels([], fontsize=x_ticklabelfontsize)

    for ax in [ax1zeta, ax2zeta, ax3zeta]:
        ax.tick_params("x", which="major", top=True, bottom=False, length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.tick_params("x", which="minor", top=True, bottom=False, length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in")
        ax.set_xlim(np.min(zeta)*0.15/np.pi, np.max(zeta)*0.15/np.pi)
        ax.set_xticks([-0.4, 0, 0.4])
        ax.set_xticks([-0.6, -0.2, 0.2, 0.6], minor=True)
    ax1zeta.set_xticklabels([r"$-2\pi/5$", r"$0$", r"$2\pi/5$"], fontsize=x_ticklabelfontsize)
    ax2zeta.set_xticklabels([], fontsize=x_ticklabelfontsize)
    ax3zeta.set_xticklabels([], fontsize=x_ticklabelfontsize)
    # for ax in [ax2, ax3]:
    #     ax.set_yticks([0, 1, 2, 3])
    #     ax.set_yticklabels([r"$0$", r"$1$", r"$2$", r"$3$"], fontsize=y_ticklabelfontsize)
    #     ax.set_yticks([0.5, 1.5, 2.5, 3.5], minor=True)
    #     ax.set_yticklabels([], minor=True)

    ax3.legend(loc="upper left", fontsize=legend_fontsize)
    # ax3.set_xticklabels([r"$0$", r"$2$", r"$4$", r"$6$", r"$8$", r"$10$"], fontsize=x_ticklabelfontsize)
    # ax1.set_xticklabels([], fontsize=x_ticklabelfontsize)
    # ax2.set_xticklabels([], fontsize=x_ticklabelfontsize)

    ax1_inset.set_xlim((-1,1))
    ax1_inset.set_ylim((-0.05, 1.05))
    ax1_inset.set_yticks([0,1])
    ax1_inset.set_yticklabels([r"$0$",r"$1$"], fontsize=inset_ticklabelfonstsize)
    ax1_inset.set_xticks([-1,0,1])
    ax1_inset.set_xticklabels([r"$-1$", r"$0$",r"$1$"], fontsize=inset_ticklabelfonstsize)
    plt.savefig("images/w003_fields_beta_scan.eps")
    plt.close()

    return

def kjm3_make_mode_plots_for_beta01_scan_for_thesis():
    """Want to scan beta for ky=1, 1.5, 2, 3 (or maybe fewer betas if time constrained)
    - For each ky, plot Omega(beta)
    - For some select ky, beta, plot:
        - phi(z)
        - apar(z)
        - bpar(z)
       together with B(z) """

    sim_longnames =  [
                make_scans.folder_name_kjm3_implicit_ky1_upw0 + "/beta_0.1000/input_for_fields",
                make_scans.folder_name_kjm3_implicit_ky1_upw0 + "/beta_0.1000/input_different_res",
                        ]
    labels = [
               r"$\beta=10\%$",
               r"$\beta=10\%$ (new resolution)",
                ]
    cols = [default_cmap(9), default_cmap(0)
            ]

    lw_list = [6, 4, 4, 2, 2]
    ls_list = ["-", "--", "-.", (0, (4,1,2,1)), "-"]

    fig = plt.figure(figsize=(12,12))
    left = 0.14
    right = 0.98
    bottom = 0.08
    top = 0.9
    vspace = 0.015
    width = right - left
    height = (top - bottom - 2*vspace)/3

    row2_bottom = bottom + height + vspace
    row1_bottom = row2_bottom + height + vspace

    x_ticklabelfontsize = 20
    x_labelfontsize = 40
    my_xticklength = 6
    my_xtickwidth = 3
    my_xticklength_minor = 3
    my_xtickwidth_minor = 1
    y_labelfontsize = x_labelfontsize
    y_ticklabelfontsize = x_ticklabelfontsize
    inset_ticklabelfonstsize = 16
    # cb_ticklabelsize = 14
    # cbax_title_fontsize = 20
    ax_title_fontsize = 20
    marker_list = ["s", "o", "P", "X", "v", "^", "<", ">", "1", "2", "3"]
    marker_size=15
    my_linewidth = 3
    legend_fontsize = 20

    inset_width = 0.26
    inset_height = 0.16
    inset_left = 0.04
    inset_bottom = 0.08

    ax1 = fig.add_axes((left, row1_bottom, width, height))
    ax2 = fig.add_axes((left, row2_bottom, width, height))
    ax3 = fig.add_axes((left, bottom, width, height))
    ax1zeta = ax1.twiny() # fig.add_axes((left, row1_bottom, width, height))
    ax2zeta = ax2.twiny() # fig.add_axes((left, row1_bottom, width, height))
    ax3zeta = ax3.twiny() # fig.add_axes((left, row1_bottom, width, height))

    for sim_idx, sim_longname in enumerate(sim_longnames):
        z, zeta, abs_phi, abs_apar, abs_bpar = get_normalised_fields_for_z_zeta(sim_longname)
        if sim_idx == 1:
            z = z * 0.2

        ax1.plot(z/np.pi, abs_phi, lw=lw_list[sim_idx], ls=ls_list[sim_idx],
                    label=labels[sim_idx], c=cols[sim_idx])
        ax2.plot(z/np.pi, abs_apar, lw=lw_list[sim_idx], ls=ls_list[sim_idx],
                    label=labels[sim_idx], c=cols[sim_idx])
        ax3.plot(z/np.pi, abs_bpar, lw=lw_list[sim_idx], ls=ls_list[sim_idx],
                    label=labels[sim_idx], c=cols[sim_idx])

    # for ax in [ax1,ax2,ax3]:
    #     ax.set_xlim((0,0.1))
    # ax1.set_ylim(0, 2.8)
    # ax2.set_ylim(-0.3, 3.9)
    # ax3.set_ylim(-0.3, 3.9)

    for ax in [ax1, ax2, ax3]:
        #ax.scatter(fprim_list, tprim_list, c="r", s=3.)
        #print("in make_plots. sorted(set(fprim_list)) = ", sorted(set(fprim_list)))
        ax.set_xlim(-0.15, 0.15)
        ax.tick_params("x", which="major", top=False, bottom=True, length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.tick_params("y", which="major", left=True, right=True, length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.tick_params("x", which="minor", top=False, bottom=True, length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in")
        ax.tick_params("y", which="minor", left=True, right=True, length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in")
        ax.set_xticks([-0.1, 0, 0.1])
        ax.set_xticks([-0.15, -0.05, 0.05, 0.15], minor=True)
        ax.set_xticklabels([], minor=True)
        ax.grid(True)

    ax3.set_xlabel(r"$z$", fontsize=x_labelfontsize)
    ax1zeta.set_xlabel(r"$\zeta$", fontsize=x_labelfontsize, labelpad=20)
    #ax1zeta.set_xlabel(r"$\zeta$", fontsize=x_labelfontsize)
    ax1.set_ylabel(r"$\vert \tilde{\phi}_k \vert$", fontsize=y_labelfontsize)
    ax2.set_ylabel(r"$\vert \tilde{A}_{\parallel k} \vert$", fontsize=y_labelfontsize)
    ax3.set_ylabel(r"$\vert \tilde{B}_{\parallel k} \vert$", fontsize=y_labelfontsize)
    ax1.set_yticks([0, 0.5, 1,])
    ax1.set_yticklabels([r"$0$", r"$0.5$", r"$1$"], fontsize=y_ticklabelfontsize)
    ax1.set_yticks([0.25, 0.75], minor=True)

    ax2.set_yticks([0, 0.05])
    ax2.set_yticklabels([r"$0$", r"$0.05$"], fontsize=y_ticklabelfontsize)
    ax2.set_yticks([0.25], minor=True)
    ax2.set_yticklabels([], minor=True)
    ax3.set_yticks([0, 0.05])
    ax3.set_yticklabels([r"$0$", r"$0.05$"], fontsize=y_ticklabelfontsize)
    ax3.set_yticks([0.25], minor=True)
    ax3.set_yticklabels([], minor=True)
    ax3.set_xticklabels([r"$-0.1$", "$0$", r"$0.1$"], fontsize=x_ticklabelfontsize)
    ax1.set_xticklabels([], fontsize=x_ticklabelfontsize)
    ax2.set_xticklabels([], fontsize=x_ticklabelfontsize)

    for ax in [ax1zeta, ax2zeta, ax3zeta]:
        ax.tick_params("x", which="major", top=True, bottom=False, length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.tick_params("x", which="minor", top=True, bottom=False, length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in")
        ax.set_xlim(np.min(zeta)*0.15/np.pi, np.max(zeta)*0.15/np.pi)
        ax.set_xticks([-0.4, 0, 0.4])
        ax.set_xticks([-0.6, -0.2, 0.2, 0.6], minor=True)
    ax1zeta.set_xticklabels([r"$-2\pi/5$", r"$0$", r"$2\pi/5$"], fontsize=x_ticklabelfontsize)
    ax2zeta.set_xticklabels([], fontsize=x_ticklabelfontsize)
    ax3zeta.set_xticklabels([], fontsize=x_ticklabelfontsize)
    # for ax in [ax2, ax3]:
    #     ax.set_yticks([0, 1, 2, 3])
    #     ax.set_yticklabels([r"$0$", r"$1$", r"$2$", r"$3$"], fontsize=y_ticklabelfontsize)
    #     ax.set_yticks([0.5, 1.5, 2.5, 3.5], minor=True)
    #     ax.set_yticklabels([], minor=True)

    ax3.legend(loc="upper left", fontsize=legend_fontsize)
    # ax3.set_xticklabels([r"$0$", r"$2$", r"$4$", r"$6$", r"$8$", r"$10$"], fontsize=x_ticklabelfontsize)
    # ax1.set_xticklabels([], fontsize=x_ticklabelfontsize)
    # ax2.set_xticklabels([], fontsize=x_ticklabelfontsize)

    plt.savefig("images/kjm3_fields_beta0.1.eps")
    plt.close()

    return

def kjm3_make_g_vpamu_colorplot_beta01_for_thesis():
    """Make a scatter plot of (vpa, mu) to show the grids"""

    def make_gvmus_plots(g_meshgrid, save_name, species):
        """ """
        my_xticklength = 5
        my_xtickwidth = 2
        x_ticklabelfontsize = 20
        y_ticklabelfontsize = 20
        y_labelfontsize = 30
        x_labelfontsize = 30
        legend_fontsize = 20
        my_linewidth = 3
        marker_size = 8
        left = 0.12
        old_right = 0.98
        top = 0.98
        bottom=0.12
        width=(old_right-left)*10/12
        right = left + width
        hspace=0.02
        height=top-bottom

        fig = plt.figure(figsize=(12,8))
        ax1 = fig.add_axes((left, bottom, width, height))

        col1_cb_lhs = right + hspace
        cb_width = 0.02
        cbax1 = fig.add_axes([col1_cb_lhs, bottom, cb_width, height])

        # fig = plt.figure(figsize=[10, 8])
        # plot_dim_x = 0.36
        # plot_dim_y = plot_dim_x*(10.0/8)
        # col1_lhs = 0.05
        # col2_lhs = 0.53
        # col2_cb_lhs = 0.93
        # row1_bottom = 0.08
        # ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
        # ax2 = fig.add_axes([col2_lhs, row1_bottom, plot_dim_x, plot_dim_y])
        #cbax2 = fig.add_axes([col2_cb_lhs, row1_bottom, cb_width, plot_dim_y])

        n_levels=100
        # if species == "ion":
        #     g_max=1
        #     g_min=
        # elif species == "electron":
        #     g_max=1
        #     g_min=0.8
        g_max=1
        g_min=np.min(g_meshgrid)
        print("g_min = ", g_min)

        g_levels = np.linspace(g_min, g_max, n_levels)
        #g_levels=n_levels
        g_contours = ax1.contourf(vpa_meshgrid, mu_meshgrid, g_meshgrid,
                                          g_levels, cmap="inferno", extend="neither")
        cbar1 = fig.colorbar(g_contours, cax=cbax1, ticklocation="left")#, ticks=[0.75, 1])
        cbax1.tick_params("y", left=False, right=True, length=my_xticklength, width=my_xtickwidth, direction="out")
        if species == "ion":
            cbar1.set_ticks([0.75,1])
            cbax1.set_yticklabels(["0.75","1"], fontsize=y_ticklabelfontsize)
        elif species == "electron":
            cbar1.set_ticks([0.8,1])
            cbax1.set_yticklabels(["0.8","1"], fontsize=y_ticklabelfontsize)

        ax1.set_ylim((np.min(mu_meshgrid), np.max(mu_meshgrid)))
        ax1.set_xlim((-3, 3))
        ax1.set_yticks([0.01, 0.1, 1])
        ax1.set_yticklabels([r"$10^{-2}$", r"$10^{-1}$", r"$10^0$"], fontsize=y_ticklabelfontsize)
        ax1.tick_params("x", length=my_xticklength, width=my_xtickwidth, direction="out")
        ax1.tick_params("y", length=my_xticklength, width=my_xtickwidth, direction="out")
        ax1.set_xticks([-3, -2, -1, 0, 1, 2, 3])
        ax1.set_xticklabels(["-3", "-2", "-1", "0", "1", "2", "3"], fontsize=x_ticklabelfontsize)
        ax1.set_ylabel(r"$\tilde{\mu}$", fontsize=y_labelfontsize)
        cbax1.yaxis.set_label_position("right")
        cbax1.yaxis.tick_right()
        if species == "ion":
            cbax1.set_ylabel(r"$\int\tilde{g}_{k,i} dz$", fontsize=y_labelfontsize, rotation=-40, labelpad=20)
        elif species == "electron":
            cbax1.set_ylabel(r"$\int\tilde{g}_{k,e} dz$", fontsize=y_labelfontsize, rotation=-40, labelpad=30)
        ax1.set_xlabel(r"$\tilde{v_\parallel}$", fontsize=x_labelfontsize)
        ax1.set_yscale("log")

        plt.savefig(save_name)
        plt.close()



        return

    outnc_longname_gvmus = "sims/kjm3_beta0.1_implicit_ky1_upw0_mode_structure/input.out.nc"
    [vpa, mu, gvmus] = extract_data_from_ncdf_with_xarray(outnc_longname_gvmus, "vpa", "mu", "gvmus")
    gvmus_tfinal_ion = np.array(gvmus[-1,0,:,:])
    gvmus_tfinal_ion = np.log(np.transpose(gvmus_tfinal_ion))
    gvmus_tfinal_electron = np.array(gvmus[-1,1,:,:])
    gvmus_tfinal_electron = np.log(np.transpose(gvmus_tfinal_electron))
    vpa_meshgrid, mu_meshgrid, [gion_meshgrid, gelectron_meshgrid] = uniquearrays2meshgrids(vpa, mu,
                [gvmus_tfinal_ion, gvmus_tfinal_electron], n_xpoints=400, n_ypoints=400)
    gion_meshgrid = gion_meshgrid/np.max(gion_meshgrid)
    gelectron_meshgrid = gelectron_meshgrid/np.max(gelectron_meshgrid)
    make_gvmus_plots(gion_meshgrid, "images/beta0.1_gvmus_vpamu_ion.png", "ion")
    make_gvmus_plots(gelectron_meshgrid, "images/beta0.1_gvmus_vpamu_electron.png", "electron")

    return

if __name__ == "__main__":
    print("Hello world")

    # make_beta_scan_plot([make_scans.folder_name_dkh_explicit_betalr,
    #                 make_scans.folder_name_dkh_explicit_betalr_fbpar0,
    #                 make_scans.folder_name_dkh_implicit_betalr,
    #                 ])
    # make_beta_scan_plot(make_scans.folder_name_dkh_explicit_betalr_fbpar0)
    # beta_scan_gene_benchmark_for_thesis()
    # make_beta_scan_plot_for_thesis()
    # make_mode_plots_for_beta_scan_for_thesis()
    # make_mode_plots_for_beta01_scan_for_thesis()
    # make_mode_plots_for_beta01_scan_for_thesis()
    # make_g_vpamu_colorplot_beta01_for_thesis()

    w003_make_beta_scan_plot_for_thesis()
    w003_make_mode_plots_for_beta_scan_for_thesis()
