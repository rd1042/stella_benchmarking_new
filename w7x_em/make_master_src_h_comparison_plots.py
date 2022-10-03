""" """

import numpy as np
import pickle
import matplotlib.pyplot as plt
import sys
sys.path.append("../postprocessing_tools")
from helper_ncdf_new import view_ncdf_variables_with_xarray, extract_data_from_ncdf_with_xarray
import plot_2d_utils as plot2dutils
from save_pickles_from_stella_scans import get_omega_data, get_phiz_data_stella

sim_longname_master_explicit = "sims/w003_es_linear_master_src_h_comparison_higher_res/master_explicit"
sim_longname_master_implicit = "sims/w003_es_linear_master_src_h_comparison_higher_res/master_str_m_implicit"
sim_longname_src_h_explicit = "sims/w003_es_linear_master_src_h_comparison_higher_res/src_h_explicit"
sim_longname_src_h_implicit = "sims/w003_es_linear_master_src_h_comparison_higher_res/src_h_str_m_implicit"
sim_longname_src_h_implicit_zupw002 = "sims/w003_es_linear_master_src_h_comparison_higher_res/src_h_str_m_implicit_zupw002"
sim_longname_src_h_implicit_tupw002 = "sims/w003_es_linear_master_src_h_comparison_higher_res/src_h_str_m_implicit_tupw002"
sim_longname_src_h_implicit_zupw002_tupw002 = "sims/w003_es_linear_master_src_h_comparison_higher_res/src_h_str_m_implicit_zupw002_tupw002"
sim_longname_master_explicit_ky35 = "sims/w003_es_linear_master_src_h_comparison_higher_res_ky3.5/master_explicit"
sim_longname_master_implicit_ky35 = "sims/w003_es_linear_master_src_h_comparison_higher_res_ky3.5/master_str_m_implicit"
sim_longname_src_h_explicit_ky35 = "sims/w003_es_linear_master_src_h_comparison_higher_res_ky3.5/src_h_explicit"
sim_longname_src_h_implicit_ky35 = "sims/w003_es_linear_master_src_h_comparison_higher_res_ky3.5/src_h_str_m_implicit"

def make_omega_t_plot_for_thesis_ky1():
    """ """
    sim_longnames = [
                    sim_longname_master_explicit,
                    sim_longname_master_implicit,
                    sim_longname_src_h_explicit,
                    # sim_longname_src_h_implicit,
                    # sim_longname_src_h_implicit_zupw002,
                    # sim_longname_src_h_implicit_tupw002,
                    sim_longname_src_h_implicit_zupw002_tupw002,
                        ]
    sim_labels = [
                 "master, explicit",
                 "master, implicit",
                  "sh, explicit",
                 "sh, implicit",
                 # "sh, implicit zup",
                 # "sh, implicit tup",
                 "sh, implicit zuptup",
    ]
    my_xticklength = 5
    my_xtickwidth = 2
    x_ticklabelfontsize = 20
    y_ticklabelfontsize = 20
    y_labelfontsize = 30
    x_labelfontsize = 30
    my_linewidth = 1#4
    left = 0.15
    right = 0.97
    top = 0.98
    bottom=0.12
    vspace = 0.02
    width=right-left
    height=(top-bottom-2*vspace)/3
    row2_bottom = bottom+height+vspace
    row1_bottom = row2_bottom+height+vspace
    ileft = 0.4
    ibottom = 0.1
    iwidth = 0.5
    iheight = 0.1

    fig = plt.figure(figsize=(10,12))

    ax1 = fig.add_axes((left, row1_bottom, width, height))
    ax2 = fig.add_axes((left, row2_bottom, width, height))
    ax3 = fig.add_axes((left, bottom, width, height))
    # inset_ax1 = fig.add_axes((ileft, row1_bottom+ibottom, iwidth, iheight))
    # inset_ax2 = fig.add_axes((ileft, row2_bottom+ibottom, iwidth, iheight))
    # inset_ax3 = fig.add_axes((ileft, bottom+ibottom, iwidth, iheight))

    for sim_idx, sim_longname in enumerate(sim_longnames):
        outnc_longname = sim_longname + ".out.nc"
        time, freqom_final, gammaom_final, freqom, gammaom, gamma_stable = get_omega_data(sim_longname, "stella")
        tend = time[-1]
        print("outnc_longname = ", outnc_longname)
        ax1.plot(time-tend, freqom, lw=my_linewidth, label=sim_labels[sim_idx])
        ax2.plot(time-tend, gammaom, lw=my_linewidth, label=sim_labels[sim_idx])
        ax3.plot(time-tend, gamma_stable, lw=my_linewidth, label=sim_labels[sim_idx])
        # inset_ax1.plot(time-tend, freqom, lw=my_linewidth)
        # inset_ax2.plot(time-tend, gammaom, lw=my_linewidth)
        # inset_ax3.plot(time-tend, gamma_stable, lw=my_linewidth)

    ax1.set_ylim((-0.4, 0.4))
    ax2.set_ylim((0, 0.272))
    ax3.set_ylim((0, 0.272))
    # inset_ax1.set_ylim((-5.5, 17))
    # inset_ax2.set_ylim((-7, 5.5))

    for ax in [ax1, ax2, ax3]:
        ax.set_xlim(-100, 0)
        ax.grid(True)
        ax.set_xticks([-100, -50, 0])
        ax.tick_params("y", length=my_xticklength, width=my_xtickwidth, direction="out")
    # for ax in [inset_ax1, inset_ax2]:
    #     ax.set_xlim(0, 200)
    #     ax.set_xticks([0, 100, 200])
    #     ax.tick_params("y", length=my_xticklength, width=my_xtickwidth, direction="in")
    #     ax.tick_params("x", length=my_xticklength, width=my_xtickwidth, direction="in")
    #     ax.set_xticklabels(["0", "100", "200"], fontsize=x_ticklabelfontsize)
    ax1.legend(loc="best")
    #ax2.legend(loc="best")
    ax1.tick_params("x", top=False, bottom=False)
    ax2.tick_params("x", top=False, bottom=True, length=my_xticklength, width=my_xtickwidth, direction="out")
    # ax1.set_yticks([0.425, 0.43, 0.435, 0.44])
    # ax1.set_yticklabels(["0.425", "0.43", "0.435", "0.44"], fontsize=y_ticklabelfontsize)
    # ax2.set_yticks([0.255, 0.26, 0.265, 0.27])
    # ax2.set_yticklabels(["0.255", "0.26", "0.265", "0.27"], fontsize=y_ticklabelfontsize)
    # inset_ax1.set_yticks([0, 10])
    # inset_ax1.set_yticklabels(["0", "10"], fontsize=y_ticklabelfontsize)
    # inset_ax2.set_yticks([-5, 0, 5])
    # inset_ax2.set_yticklabels(["-5", "0", "5"], fontsize=y_ticklabelfontsize)
    # ax1.set_xticklabels([])
    # ax2.set_xticklabels(["100", "150", "200"], fontsize=x_ticklabelfontsize)
    ax1.set_ylabel(r"$\tilde{\omega}$", fontsize=y_labelfontsize)
    ax2.set_ylabel(r"$\tilde{\gamma}$", fontsize=y_labelfontsize)
    ax2.set_xlabel(r"$\tilde{t}$", fontsize=x_labelfontsize)
    # plt.show()
    plt.savefig("images/master_src_h_omega_t.eps")
    plt.close()
    return

def make_omega_t_plot_for_thesis_ky35():
    """ """
    sim_longnames = [
                sim_longname_master_explicit_ky35,
                sim_longname_master_implicit_ky35,
                sim_longname_src_h_explicit_ky35,
                sim_longname_src_h_implicit_ky35
                        ]
    sim_labels = [
                 "master, explicit",
                 "master, implicit",
                  "sh, explicit",
                 "sh, implicit",
                ]
    my_xticklength = 5
    my_xtickwidth = 2
    x_ticklabelfontsize = 20
    y_ticklabelfontsize = 20
    y_labelfontsize = 30
    x_labelfontsize = 30
    my_linewidth = 1#4
    left = 0.15
    right = 0.97
    top = 0.98
    bottom=0.12
    vspace = 0.02
    width=right-left
    height=(top-bottom-2*vspace)/3
    row2_bottom = bottom+height+vspace
    row1_bottom = row2_bottom+height+vspace
    ileft = 0.4
    ibottom = 0.1
    iwidth = 0.5
    iheight = 0.1

    fig = plt.figure(figsize=(10,12))

    ax1 = fig.add_axes((left, row1_bottom, width, height))
    ax2 = fig.add_axes((left, row2_bottom, width, height))
    ax3 = fig.add_axes((left, bottom, width, height))
    # inset_ax1 = fig.add_axes((ileft, row1_bottom+ibottom, iwidth, iheight))
    # inset_ax2 = fig.add_axes((ileft, row2_bottom+ibottom, iwidth, iheight))
    # inset_ax3 = fig.add_axes((ileft, bottom+ibottom, iwidth, iheight))

    for sim_idx, sim_longname in enumerate(sim_longnames):
        outnc_longname = sim_longname + ".out.nc"
        time, freqom_final, gammaom_final, freqom, gammaom, gamma_stable = get_omega_data(sim_longname, "stella")
        tend = time[-1]
        print("outnc_longname = ", outnc_longname)
        ax1.plot(time-tend, freqom, lw=my_linewidth, label=sim_labels[sim_idx])
        ax2.plot(time-tend, gammaom, lw=my_linewidth, label=sim_labels[sim_idx])
        ax3.plot(time-tend, gamma_stable, lw=my_linewidth, label=sim_labels[sim_idx])
        # inset_ax1.plot(time-tend, freqom, lw=my_linewidth)
        # inset_ax2.plot(time-tend, gammaom, lw=my_linewidth)
        # inset_ax3.plot(time-tend, gamma_stable, lw=my_linewidth)

    ax1.set_ylim((-0.4, 0.4))
    ax2.set_ylim((0, 0.272))
    ax3.set_ylim((0, 0.272))
    # inset_ax1.set_ylim((-5.5, 17))
    # inset_ax2.set_ylim((-7, 5.5))

    for ax in [ax1, ax2, ax3]:
        ax.set_xlim(-100, 0)
        ax.grid(True)
        ax.set_xticks([-100, -50, 0])
        ax.tick_params("y", length=my_xticklength, width=my_xtickwidth, direction="out")
    # for ax in [inset_ax1, inset_ax2]:
    #     ax.set_xlim(0, 200)
    #     ax.set_xticks([0, 100, 200])
    #     ax.tick_params("y", length=my_xticklength, width=my_xtickwidth, direction="in")
    #     ax.tick_params("x", length=my_xticklength, width=my_xtickwidth, direction="in")
    #     ax.set_xticklabels(["0", "100", "200"], fontsize=x_ticklabelfontsize)
    ax1.legend(loc="best")
    #ax2.legend(loc="best")
    ax1.tick_params("x", top=False, bottom=False)
    ax2.tick_params("x", top=False, bottom=True, length=my_xticklength, width=my_xtickwidth, direction="out")
    # ax1.set_yticks([0.425, 0.43, 0.435, 0.44])
    # ax1.set_yticklabels(["0.425", "0.43", "0.435", "0.44"], fontsize=y_ticklabelfontsize)
    # ax2.set_yticks([0.255, 0.26, 0.265, 0.27])
    # ax2.set_yticklabels(["0.255", "0.26", "0.265", "0.27"], fontsize=y_ticklabelfontsize)
    # inset_ax1.set_yticks([0, 10])
    # inset_ax1.set_yticklabels(["0", "10"], fontsize=y_ticklabelfontsize)
    # inset_ax2.set_yticks([-5, 0, 5])
    # inset_ax2.set_yticklabels(["-5", "0", "5"], fontsize=y_ticklabelfontsize)
    # ax1.set_xticklabels([])
    # ax2.set_xticklabels(["100", "150", "200"], fontsize=x_ticklabelfontsize)
    ax1.set_ylabel(r"$\tilde{\omega}$", fontsize=y_labelfontsize)
    ax2.set_ylabel(r"$\tilde{\gamma}$", fontsize=y_labelfontsize)
    ax2.set_xlabel(r"$\tilde{t}$", fontsize=x_labelfontsize)
    # plt.show()
    plt.savefig("images/master_src_h_omega_t.eps")
    plt.close()
    return

def make_phi_z_plot_for_thesis_ky35():
    """ """

    my_xticklength = 5
    my_xtickwidth = 2
    x_ticklabelfontsize = 20
    y_ticklabelfontsize = 20
    y_labelfontsize = 30
    x_labelfontsize = 30
    legend_fontsize = 20
    my_linewidth = 3
    left = 0.2
    right = 0.9
    top = 0.86
    bottom=0.12
    width=right-left
    height=top-bottom

    fig = plt.figure(figsize=(12,8))
    ax3 = fig.add_axes((left, bottom, width, height))
    ax2 = ax3.twiny()
    ax1 = ax3.twinx()

    sim_longnames = [
                sim_longname_master_explicit_ky35,
                sim_longname_master_implicit_ky35,
                sim_longname_src_h_explicit_ky35,
                sim_longname_src_h_implicit_ky35
                        ]
    sim_labels = [
                 "master, explicit",
                 "master, implicit",
                  "sh, explicit",
                 "sh, implicit",
                ]

    for sim_idx, sim_longname in enumerate(sim_longnames):
        z, real_phi, imag_phi = get_phiz_data_stella(sim_longname)
        # z = np.array(z)

        geo_longname = sim_longname + ".geometry"
        data = np.loadtxt(geo_longname, skiprows=4, usecols=(1,2,3)) # data = [zed, zeta, bmaf]
        zed = data[:,0]
        zeta = data[:,1]
        bmag = data[:,2]
        min_zeta = np.min(zeta)
        max_zeta = np.max(zeta)
        abs_phi = np.sqrt(real_phi*real_phi + imag_phi*imag_phi)
        max_abs_phi = np.max(abs_phi)
        real_phi = real_phi/max_abs_phi
        imag_phi = imag_phi/max_abs_phi
        abs_phi = abs_phi/max_abs_phi

        if sim_idx == 1:
            ax3.plot(zed/np.pi, bmag, label = r"$\tilde{B}$",c="gray", ls="--", lw=1, alpha=0.5, zorder=10)
        # ax1.plot(z/np.pi, real_phi, label=r"re$(\tilde{\varphi}_k)$", lw=my_linewidth, zorder=100)
        # ax1.plot(z/np.pi, imag_phi, label=r"im$(\tilde{\varphi}_k)$", lw=my_linewidth, zorder=100)
        ax1.plot(z/np.pi, abs_phi, label=sim_labels[sim_idx], lw=my_linewidth, zorder=100)

    ax3.yaxis.set_label_position("right")
    ax1.yaxis.set_label_position("left")
    ax3.yaxis.tick_right()
    ax1.yaxis.tick_left()
    ax1.plot([-10, -11], [-10, -11], label = r"$\tilde{B}$",c="gray", ls="--", lw=1, alpha=0.5)
    ax1.legend(loc="best", fontsize=legend_fontsize)
    ax1.set_ylim((-0.1, 1.1))
    ax1.set_xlim((-1, 1))
    ax3.set_ylim((0.95, 1.22))
    ax2.set_xlim((min_zeta/np.pi, max_zeta/np.pi))
    ax1.set_yticks([-0.5, 0, 0.5, 1])
    ax1.set_yticklabels([r"$-0.5$", "$0$", r"$0.5$", r"$1$"], fontsize=y_ticklabelfontsize)
    ax3.set_yticks([1, 1.1, 1.2])
    ax3.set_yticklabels([r"$1$", "$1.1$", r"$1.2$"], fontsize=y_ticklabelfontsize)
    ax3.tick_params("x", top=False, bottom=True, length=my_xticklength, width=my_xtickwidth, direction="out")
    ax3.set_xticks([-1, 0, 1])
    ax3.set_xticklabels([r"$-\pi$", "0", r"$\pi$"], fontsize=x_ticklabelfontsize)
    ax2.tick_params("x", top=True, bottom=False, length=my_xticklength, width=my_xtickwidth, direction="out")
    ax2.set_xticks([-4,  0, 4])
    ax2.set_xticklabels([r"$-4\pi$", "0", r"$4\pi$"], fontsize=x_ticklabelfontsize)
    ax1.set_ylabel(r"$\frac{\tilde{\varphi}_k}{max(\vert \tilde{\varphi}_k \vert)}$", fontsize=y_labelfontsize,
                    labelpad=20, rotation=30)
    ax3.set_ylabel(r"$\tilde{B}$", fontsize=y_labelfontsize,
                    labelpad=10, rotation=0)
    ax1.set_xlabel(r"$z$", fontsize=x_labelfontsize)
    ax2.set_xlabel(r"$\zeta$", fontsize=x_labelfontsize, labelpad=20)


    #ax1.grid(True)
    plt.savefig("images/master_src_h_comparison_phi_z.eps")
    plt.close()
    return

if __name__ == "__main__":
    print("Hello world")
    make_omega_t_plot_for_thesis_ky35()
    make_phi_z_plot_for_thesis_ky35()
