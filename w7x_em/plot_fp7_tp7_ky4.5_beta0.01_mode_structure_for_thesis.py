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

def w003_make_mode_plots_for_sim():
    """Want to scan beta for ky=1, 1.5, 2, 3 (or maybe fewer betas if time constrained)
    - For each ky, plot Omega(beta)
    - For some select ky, beta, plot:
        - phi(z)
        - apar(z)
        - bpar(z)
       together with B(z) """

    sim_longnames =  [
                    "sims/w003_fprim_7.7143_tprim_7.7143_ky_4.5000_beta0.01_for_diagnostics/input_short_time"
                        ]
    labels = [ "",
                ""
                ]
    cols = [default_cmap(0), default_cmap(7),
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

    ax2.set_yticks([0, 0.1, 0.2])
    ax2.set_yticklabels([r"$0$", r"$0.1$", r"$0.2$"], fontsize=y_ticklabelfontsize)
    ax2.set_yticks([0.05, 0.15], minor=True)
    ax2.set_yticklabels([], minor=True)
    ax3.set_yticks([0, 0.01, 0.02])
    ax3.set_yticklabels([r"$0$", r"$0.01$", r"$0.02$"], fontsize=y_ticklabelfontsize)
    ax3.set_yticks([0.015, 0.005], minor=True)
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

    # ax3.legend(loc="upper left", fontsize=legend_fontsize)
    # ax3.set_xticklabels([r"$0$", r"$2$", r"$4$", r"$6$", r"$8$", r"$10$"], fontsize=x_ticklabelfontsize)
    # ax1.set_xticklabels([], fontsize=x_ticklabelfontsize)
    # ax2.set_xticklabels([], fontsize=x_ticklabelfontsize)

    ax1_inset.set_xlim((-1,1))
    ax1_inset.set_ylim((-0.05, 1.05))
    ax1_inset.set_yticks([0,1])
    ax1_inset.set_yticklabels([r"$0$",r"$1$"], fontsize=inset_ticklabelfonstsize)
    ax1_inset.set_xticks([-1,0,1])
    ax1_inset.set_xticklabels([r"$-1$", r"$0$",r"$1$"], fontsize=inset_ticklabelfonstsize)
    plt.savefig("images/w003_mode_plot_fprim77_tprim77_ky45.eps")
    plt.close()

    return

def w003_make_mode_plots_for_sim_plus_extras():
    """Want to scan beta for ky=1, 1.5, 2, 3 (or maybe fewer betas if time constrained)
    - For each ky, plot Omega(beta)
    - For some select ky, beta, plot:
        - phi(z)
        - apar(z)
        - bpar(z)
       together with B(z) """

    sim_longnames =  [
                    "sims/w003_fprim_7.7143_tprim_7.7143_ky_4.5000_beta0.01_for_diagnostics/input_short_time",
                    "sims/w003_fprim_7.7143_tprim_7.7143_ky_4.5000_beta0.01_for_diagnostics/input_short_time_all_fprim",
                    "sims/w003_fprim_7.7143_tprim_7.7143_ky_4.5000_beta0.01_for_diagnostics/input_short_time_all_tprim",
                    "sims/w003_fprim_7.7143_tprim_7.7143_ky_4.5000_beta0.01_for_diagnostics/input_short_time_beta0.02",
                        ]
    labels = ["split",
              "all fprim",
              "all tprim",
              "beta=2%",
                ]
    cols = [default_cmap(0), default_cmap(7),
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
    ax2.legend(loc="best")
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

    plt.show()
    sys.exit()
    ax3.set_xlabel(r"$z$", fontsize=x_labelfontsize)
    ax1zeta.set_xlabel(r"$\zeta$", fontsize=x_labelfontsize, labelpad=20)
    #ax1zeta.set_xlabel(r"$\zeta$", fontsize=x_labelfontsize)
    ax1.set_ylabel(r"$\vert \tilde{\phi}_k \vert$", fontsize=y_labelfontsize)
    ax2.set_ylabel(r"$\vert \tilde{A}_{\parallel k} \vert$", fontsize=y_labelfontsize)
    ax3.set_ylabel(r"$\vert \tilde{B}_{\parallel k} \vert$", fontsize=y_labelfontsize)
    ax1.set_yticks([0, 0.5, 1,])
    ax1.set_yticklabels([r"$0$", r"$0.5$", r"$1$"], fontsize=y_ticklabelfontsize)
    ax1.set_yticks([0.25, 0.75], minor=True)

    ax2.set_yticks([0, 0.1, 0.2])
    ax2.set_yticklabels([r"$0$", r"$0.1$", r"$0.2$"], fontsize=y_ticklabelfontsize)
    ax2.set_yticks([0.05, 0.15], minor=True)
    ax2.set_yticklabels([], minor=True)
    ax3.set_yticks([0, 0.01, 0.02])
    ax3.set_yticklabels([r"$0$", r"$0.01$", r"$0.02$"], fontsize=y_ticklabelfontsize)
    ax3.set_yticks([0.015, 0.005], minor=True)
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

    # ax3.legend(loc="upper left", fontsize=legend_fontsize)
    # ax3.set_xticklabels([r"$0$", r"$2$", r"$4$", r"$6$", r"$8$", r"$10$"], fontsize=x_ticklabelfontsize)
    # ax1.set_xticklabels([], fontsize=x_ticklabelfontsize)
    # ax2.set_xticklabels([], fontsize=x_ticklabelfontsize)

    ax1_inset.set_xlim((-1,1))
    ax1_inset.set_ylim((-0.05, 1.05))
    ax1_inset.set_yticks([0,1])
    ax1_inset.set_yticklabels([r"$0$",r"$1$"], fontsize=inset_ticklabelfonstsize)
    ax1_inset.set_xticks([-1,0,1])
    ax1_inset.set_xticklabels([r"$-1$", r"$0$",r"$1$"], fontsize=inset_ticklabelfonstsize)
    plt.savefig("images/w003_mode_plot_fprim77_tprim77_ky45_plus_scans.eps")
    plt.close()

    return

if __name__ == "__main__":
    print("Hello world")
    # w003_make_mode_plots_for_sim()
    w003_make_mode_plots_for_sim_plus_extras()
