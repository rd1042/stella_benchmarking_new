""" """

import sys
sys.path.append("../postprocessing_tools")
import matplotlib.scale as scale
import matplotlib.pyplot as plt
import numpy as np
import glob
import re
import pickle
from scipy.interpolate import griddata

IMAGE_DIR = "./images/"

master_pickle = "sims/stella_nonlinear_2species_master_archer2/input.summary_pickle"
master_cfl025_pickle = "sims/stella_nonlinear_2species_master_cflcushion025/input.summary_pickle"
leapfrog_pickle = "sims/stella_nonlinear_2species_leapfrog_nonlinear/input.summary_pickle"
nisl_pickle = "sims/stella_nonlinear_2species_nisl_archer2/input.summary_pickle"
isl_pickle = "sims/stella_nonlinear_2species_isl/input.summary_pickle"
default_cmap = plt.get_cmap("tab10")

def plot_phi2_t(pickle_longnames):
    """ """

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    for pickle_longname in pickle_longnames:
        shortname = re.split("/", pickle_longname)[-2]
        myfile = open(pickle_longname, "rb")
        [t, kx, ky, phi2, phi2_kxky] = pickle.load(myfile)
        myfile.close()
        ax1.plot(t, phi2, label=shortname, lw=2)

    ax1.set_yscale("log")
    ax1.set_xlabel(r"$t$")
    ax1.set_ylabel(r"phi2")
    ax1.legend(loc="best")
    plt.show()

    return

def make_colorplot(fig, ax, cbax, kx, ky, phi2_kxky, zonal=True, log=True):
    """ """
    if zonal:
        phi2 = (np.swapaxes(phi2_kxky, 0, 1))

    else:
        phi2 = (np.swapaxes(phi2_kxky[:,1:], 0, 1))

    print("scale names:", scale.get_scale_names())
    if log:
        img = ax.imshow(np.log(phi2), aspect="auto", origin="lower")
    else:
        img = ax.imshow(phi2, aspect="auto", origin="lower")
    cbar = fig.colorbar(img, cax=cbax)

    return cbar

def plot_properties_master_sim(cfl_cushion="0.5"):
    """ """
    if cfl_cushion == "0.5":
        pickle = master_pickle
    else:
        pickle = master_cfl025_pickle

    [t, kx, ky, phi2, phi2_tkxky, parsed_pflx, parsed_vflx, parsed_qflx] = get_data_from_pickle(pickle, fluxes=True)

    ## Dimensions of phi2_tkxky are [t, kx, ky]
    ## Dimensions of parsed_pflx are [t_idxs, spec, kx, ky]

    t_idxs_sample = [0, 11, 100, -1]

    ## imshow makes the plot, but in index-space rather than (kx, ky)
    # To have meaningful labels, we want to find indexes corresponding to the
    # desired values of kx, ky
    kx_label_vals = np.array([-2, -1, 0, 1, 2])
    ky_label_vals = np.array([0, 1, 2])
    # idx space goes [0, 1, 2, . . . , nkx(y)-1]
    # kx(y) space goes [kx(y)_min,  . . . , kx(y)_max]
    # so kx = ((kx_max - kx_min)/(nkx-1))*kxidx + kx_min
    # so kxidx = (kx - kx_min) * ((nkx-1)/(kx_max - kx_min))
    kx_min = np.min(kx) ; kx_max = np.max(kx)
    ky_min = np.min(ky) ; ky_max = np.max(ky)

    kx_extent = kx_max - kx_min
    ky_extent = ky_max - ky_min
    nkx = len(kx) ; nky = len(ky)
    xlabel_idx_vals = (kx_label_vals - kx_min) * ((nkx-1)/kx_extent)
    ylabel_idx_vals = (ky_label_vals - ky_min) * ((nky-1)/ky_extent)
    # print("xlabel_idx_vals = ", xlabel_idx_vals)
    # print("ylabel_idx_vals = ", ylabel_idx_vals)
    # sys.exit()
    marker_size = 50
    marker_size_large = 300
    linewidth=2

    left = 0.16
    right =0.95
    top = 0.95
    bottom = 0.14
    width = right - left
    height = top - bottom
    fig = plt.figure(figsize=(12,10))
    ax1 = fig.add_axes((left, bottom, width, height))

    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    ax1.plot(t, phi2, lw=linewidth, zorder=0, c="black", )
    ax1.scatter(t, phi2, c="red", s=marker_size, zorder=10, marker="x")
    for t_idx in t_idxs_sample:
        ax1.scatter(t[t_idx], phi2[t_idx], c="blue", s=marker_size_large, zorder=10, marker="X")

    ax1.set_yscale("log")

    legend_fontsize = 16
    xlabel_fontsize = 30
    ylabel_fontsize = 30
    xticklabel_fontsize = 18
    yticklabel_fontsize = 18
    ax1.set_xlabel(r"$\tilde{t}$", fontsize=xlabel_fontsize)
    ax1.set_ylabel(r"$ \sum_{\tilde{k}_x, \tilde{k}_y} \vert \tilde{\phi}_k(t) \vert^2$", fontsize=ylabel_fontsize)
    ax1.set_xlim(0, 450)
    ax1.set_ylim(2E-4, 420)
    ax1.set_xticks([0, 100, 200, 300, 400])
    ax1.set_xticklabels([r"$0$", r"$100$", r"$200$", r"$300$", r"$400$"], fontsize=xticklabel_fontsize)
    ax1.set_yticks([0.01, 1, 100,])
    ax1.set_yticklabels([r"$10^{-2}$", r"$10^0$", r"$10^2$"], fontsize=yticklabel_fontsize)



    if cfl_cushion == "0.5":
        plt.savefig("master_phi2_t.eps")
    else:
        plt.savefig("master_phi2_t_cushion025.eps")
    plt.close()

    fig = plt.figure(figsize=[12, 12])

    left = 0.1
    right = 0.95
    top = 0.92
    bottom = 0.1

    # row1_bottom = 0.1
    subplot_hspace = 0.1
    subplot_vspace = 0.1
    subplot_width = (right - left - subplot_hspace)/2
    subplot_height = (top - bottom - subplot_vspace)/2
    cb_width = 0.02
    cb_hspace = 0.02
    plot_width = subplot_width - cb_width - cb_hspace
    plot_height = subplot_height
    print("plot_width, plot_height = ", plot_width, plot_height)

    ax_col1_lhs = left
    ax_col2_lhs = left + subplot_width + subplot_hspace
    cbax_col1_lhs = left + plot_width + cb_hspace
    cbax_col2_lhs = left + subplot_width + subplot_hspace + plot_width + cb_hspace
    row1_bottom = bottom + subplot_vspace + subplot_height

    xticklabel_fontsize = 24
    yticklabel_fontsize = 24
    xlabel_fontsize = 30
    ylabel_fontsize = 30
    cbax_title_fontsize = 30
    # plot_dim_x = 0.75
    # plot_dim_y = plot_dim_x
    # col1_cb_lhs = 0.9

    ax1 = fig.add_axes([ax_col1_lhs, row1_bottom, plot_width, plot_height])
    cbax1 = fig.add_axes([cbax_col1_lhs, row1_bottom, cb_width, plot_height])
    ax2 = fig.add_axes([ax_col2_lhs, row1_bottom, plot_width, plot_height])
    cbax2 = fig.add_axes([cbax_col2_lhs, row1_bottom, cb_width, plot_height])
    ax3 = fig.add_axes([ax_col1_lhs, bottom, plot_width, plot_height])
    cbax3 = fig.add_axes([cbax_col1_lhs, bottom, cb_width, plot_height])
    ax4 = fig.add_axes([ax_col2_lhs, bottom, plot_width, plot_height])
    cbax4 = fig.add_axes([cbax_col2_lhs, bottom, cb_width, plot_height])

    axes = [ax1, ax2, ax3, ax4]
    cbaxes = [cbax1, cbax2, cbax3, cbax4]
    cbar_list = []
    for counter, t_idx in enumerate(t_idxs_sample):
        phi2_kxky = phi2_tkxky[t_idx,:,:]
        ax = axes[counter]  ; cbax = cbaxes[counter]
        cbar = make_colorplot(fig, ax, cbax, kx, ky, np.log10(phi2_kxky), zonal=True, log=False)
        cbar_list.append(cbar)

    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_xticks(xlabel_idx_vals)
        ax.set_yticks(ylabel_idx_vals)
        ax.set_xticklabels([r"$-2$", r"$-1$", r"$0$", r"$1$", r"$2$", ], fontsize=xticklabel_fontsize)
        ax.set_yticklabels([r"$0$", r"$1$", r"$2$"], fontsize=yticklabel_fontsize)

    cbarticklabel_fontsize = 16
    cbar1 = cbar_list[0]
    cbar2 = cbar_list[1]
    cbar3 = cbar_list[2]
    cbar4 = cbar_list[3]
    cbar1.set_ticks([-250, -200, -150, -100, -50])
    cbax1.set_yticklabels([r"$-250$", r"$-200$", r"$-150$", r"$-100$", r"$-50$"], fontsize=cbarticklabel_fontsize)
    cbar2.set_ticks([-8, -6, -4])
    cbax2.set_yticklabels([r"$-8$", r"$-6$", r"$-4$"], fontsize=cbarticklabel_fontsize)
    cbar3.set_ticks([-3, -2, -1, 0, 1])
    cbax3.set_yticklabels([r"$-3$", r"$-2$", r"$-1$", r"$0$", r"$1$"], fontsize=cbarticklabel_fontsize)
    cbar4.set_ticks([-3, -2, -1, 0])
    cbax4.set_yticklabels([r"$-3$", r"$-2$", r"$-1$", r"$0$"], fontsize=cbarticklabel_fontsize)
    for ax in [ax3, ax4]:
        ax.set_xlabel(r"$\tilde{k}_x$", fontsize=xlabel_fontsize)
    for ax in [ax1, ax3]:
        ax.set_ylabel(r"$\tilde{k}_y$", fontsize=ylabel_fontsize)

    for cbax in [cbax1, cbax2]:
        cbax.set_title(r"$\vert \tilde{\phi}_k\vert^2$", fontsize=cbax_title_fontsize)

    if cfl_cushion == "0.5":
        plt.savefig("phi2_map_master.eps")
    else:
        plt.savefig("phi2_map_master_cushion025.eps")
    plt.close()

    ########################################################################
    ##################### FLUXES ###########################################

    fig = plt.figure(figsize=[12, 12])
    ax1 = fig.add_axes([ax_col1_lhs, row1_bottom, plot_width, plot_height])
    cbax1 = fig.add_axes([cbax_col1_lhs, row1_bottom, cb_width, plot_height])
    ax2 = fig.add_axes([ax_col2_lhs, row1_bottom, plot_width, plot_height])
    cbax2 = fig.add_axes([cbax_col2_lhs, row1_bottom, cb_width, plot_height])
    ax3 = fig.add_axes([ax_col1_lhs, bottom, plot_width, plot_height])
    cbax3 = fig.add_axes([cbax_col1_lhs, bottom, cb_width, plot_height])
    ax4 = fig.add_axes([ax_col2_lhs, bottom, plot_width, plot_height])
    cbax4 = fig.add_axes([cbax_col2_lhs, bottom, cb_width, plot_height])

    axes = [ax1, ax2, ax3, ax4]
    cbaxes = [cbax1, cbax2, cbax3, cbax4]
    cbar = make_colorplot(fig, ax1, cbax1, kx, ky, parsed_qflx[-1,0,:,:], zonal=True, log=False)
    cbar = make_colorplot(fig, ax2, cbax2, kx, ky, parsed_qflx[-1,1,:,:], zonal=True, log=False)
    cbar = make_colorplot(fig, ax3, cbax3, kx, ky, parsed_pflx[-1,0,:,:], zonal=True, log=False)
    cbar = make_colorplot(fig, ax4, cbax4, kx, ky, parsed_pflx[-1,1,:,:], zonal=True, log=False)

    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_xticks(xlabel_idx_vals)
        ax.set_yticks(ylabel_idx_vals)
        ax.set_xticklabels([r"$-2$", r"$-1$", r"$0$", r"$1$", r"$2$", ], fontsize=xticklabel_fontsize)
        ax.set_yticklabels([r"$0$", r"$1$", r"$2$"], fontsize=yticklabel_fontsize)
        ax.set_xlabel(r"$k_x$", fontsize=xlabel_fontsize)
        ax.set_ylabel(r"$k_y$", fontsize=ylabel_fontsize)

    for cbax in [cbax1, cbax2]:
        cbax.set_title(r"$\vert \tilde{\phi}_k\vert^2$", fontsize=cbax_title_fontsize)

    if cfl_cushion == "0.5":
        plt.savefig("fluxes_master.eps")
    else:
        plt.savefig("fluxes_master_cushion025.eps")
    plt.close()

    return

def make_master_cushion_comparison():
    """ """
    [t_orig, kx, ky, phi2_orig, phi2_tkxky_orig, parsed_pflx_orig, parsed_vflx_orig, parsed_qflx_orig] = get_data_from_pickle(master_pickle, fluxes=True)
    [t_cfl025, kx, ky, phi2_cfl025, phi2_tkxky_cfl025, parsed_pflx_cfl025, parsed_vflx_cfl025, parsed_qflx_cfl025] = get_data_from_pickle(master_cfl025_pickle, fluxes=True)

    ## Dimensions of phi2_tkxky are [t, kx, ky]
    ## Dimensions of parsed_pflx are [t_idxs, spec, kx, ky]

    t_idxs_sample = [0, 11, 100, -1]

    ## imshow makes the plot, but in index-space rather than (kx, ky)
    # To have meaningful labels, we want to find indexes corresponding to the
    # desired values of kx, ky
    kx_label_vals = np.array([-2, -1, 0, 1, 2])
    ky_label_vals = np.array([0, 1, 2])
    # idx space goes [0, 1, 2, . . . , nkx(y)-1]
    # kx(y) space goes [kx(y)_min,  . . . , kx(y)_max]
    # so kx = ((kx_max - kx_min)/(nkx-1))*kxidx + kx_min
    # so kxidx = (kx - kx_min) * ((nkx-1)/(kx_max - kx_min))
    kx_min = np.min(kx) ; kx_max = np.max(kx)
    ky_min = np.min(ky) ; ky_max = np.max(ky)

    kx_extent = kx_max - kx_min
    ky_extent = ky_max - ky_min
    nkx = len(kx) ; nky = len(ky)
    xlabel_idx_vals = (kx_label_vals - kx_min) * ((nkx-1)/kx_extent)
    ylabel_idx_vals = (ky_label_vals - ky_min) * ((nky-1)/ky_extent)
    # print("xlabel_idx_vals = ", xlabel_idx_vals)
    # print("ylabel_idx_vals = ", ylabel_idx_vals)
    # sys.exit()
    marker_size = 5
    linewidth=2
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(t_orig, phi2_orig, lw=linewidth, zorder=0)
    ax1.scatter(t_orig, phi2_orig, c="red", s=marker_size, zorder=10)
    ax1.plot(t_cfl025, phi2_cfl025, lw=linewidth, zorder=0)
    ax1.scatter(t_cfl025, phi2_cfl025, c="red", s=marker_size, zorder=10)
    for t_idx in t_idxs_sample:
        ax1.scatter(t_orig[t_idx], phi2_orig[t_idx], c="blue", s=marker_size, zorder=10)
        ax1.scatter(t_cfl025[t_idx], phi2_cfl025[t_idx], c="blue", s=marker_size, zorder=10)

    ax1.set_yscale("log")
    ax1.set_xlabel(r"$t$")
    ax1.set_ylabel(r"phi2")
    plt.tight_layout()
    plt.savefig("master_dtscan_phi2_t.eps")
    plt.close()

    # fig = plt.figure(figsize=[12, 12])
    #
    # left = 0.1
    # right = 0.95
    # top = 0.92
    # bottom = 0.1
    #
    # # row1_bottom = 0.1
    # subplot_hspace = 0.1
    # subplot_vspace = 0.1
    # subplot_width = (right - left - subplot_hspace)/2
    # subplot_height = (top - bottom - subplot_vspace)/2
    # cb_width = 0.02
    # cb_hspace = 0.02
    # plot_width = subplot_width - cb_width - cb_hspace
    # plot_height = subplot_height
    # print("plot_width, plot_height = ", plot_width, plot_height)
    #
    # ax_col1_lhs = left
    # ax_col2_lhs = left + subplot_width + subplot_hspace
    # cbax_col1_lhs = left + plot_width + cb_hspace
    # cbax_col2_lhs = left + subplot_width + subplot_hspace + plot_width + cb_hspace
    # row1_bottom = bottom + subplot_vspace + subplot_height
    #
    # xticklabel_fontsize = 24
    # yticklabel_fontsize = 24
    # xlabel_fontsize = 30
    # ylabel_fontsize = 30
    # cbax_title_fontsize = 30
    # # plot_dim_x = 0.75
    # # plot_dim_y = plot_dim_x
    # # col1_cb_lhs = 0.9
    #
    # ax1 = fig.add_axes([ax_col1_lhs, row1_bottom, plot_width, plot_height])
    # cbax1 = fig.add_axes([cbax_col1_lhs, row1_bottom, cb_width, plot_height])
    # ax2 = fig.add_axes([ax_col2_lhs, row1_bottom, plot_width, plot_height])
    # cbax2 = fig.add_axes([cbax_col2_lhs, row1_bottom, cb_width, plot_height])
    # ax3 = fig.add_axes([ax_col1_lhs, bottom, plot_width, plot_height])
    # cbax3 = fig.add_axes([cbax_col1_lhs, bottom, cb_width, plot_height])
    # ax4 = fig.add_axes([ax_col2_lhs, bottom, plot_width, plot_height])
    # cbax4 = fig.add_axes([cbax_col2_lhs, bottom, cb_width, plot_height])
    #
    # axes = [ax1, ax2, ax3, ax4]
    # cbaxes = [cbax1, cbax2, cbax3, cbax4]
    # for counter, t_idx in enumerate(t_idxs_sample):
    #     phi2_kxky = phi2_tkxky[t_idx,:,:]
    #     ax = axes[counter]  ; cbax = cbaxes[counter]
    #     make_colorplot(fig, ax, cbax, kx, ky, phi2_kxky, zonal=True)
    #
    # for ax in [ax1, ax2, ax3, ax4]:
    #     ax.set_xticks(xlabel_idx_vals)
    #     ax.set_yticks(ylabel_idx_vals)
    #     ax.set_xticklabels([r"$-2$", r"$-1$", r"$0$", r"$1$", r"$2$", ], fontsize=xticklabel_fontsize)
    #     ax.set_yticklabels([r"$0$", r"$1$", r"$2$"], fontsize=yticklabel_fontsize)
    #     ax.set_xlabel(r"$k_x$", fontsize=xlabel_fontsize)
    #     ax.set_ylabel(r"$k_y$", fontsize=ylabel_fontsize)
    #
    # for cbax in [cbax1, cbax2]:
    #     cbax.set_title(r"$\vert \tilde{\phi}_k\vert^2$", fontsize=cbax_title_fontsize)
    #
    # plt.savefig("phi2_map_master.eps")
    # plt.close()
    #
    #
    # fig = plt.figure(figsize=[12, 12])
    # ax1 = fig.add_axes([ax_col1_lhs, row1_bottom, plot_width, plot_height])
    # cbax1 = fig.add_axes([cbax_col1_lhs, row1_bottom, cb_width, plot_height])
    # ax2 = fig.add_axes([ax_col2_lhs, row1_bottom, plot_width, plot_height])
    # cbax2 = fig.add_axes([cbax_col2_lhs, row1_bottom, cb_width, plot_height])
    # ax3 = fig.add_axes([ax_col1_lhs, bottom, plot_width, plot_height])
    # cbax3 = fig.add_axes([cbax_col1_lhs, bottom, cb_width, plot_height])
    # ax4 = fig.add_axes([ax_col2_lhs, bottom, plot_width, plot_height])
    # cbax4 = fig.add_axes([cbax_col2_lhs, bottom, cb_width, plot_height])
    #
    # axes = [ax1, ax2, ax3, ax4]
    # cbaxes = [cbax1, cbax2, cbax3, cbax4]
    # make_colorplot(fig, ax1, cbax1, kx, ky, parsed_qflx[-1,0,:,:], zonal=True, log=False)
    # make_colorplot(fig, ax2, cbax2, kx, ky, parsed_qflx[-1,1,:,:], zonal=True, log=False)
    # make_colorplot(fig, ax3, cbax3, kx, ky, parsed_pflx[-1,0,:,:], zonal=True, log=False)
    # make_colorplot(fig, ax4, cbax4, kx, ky, parsed_pflx[-1,1,:,:], zonal=True, log=False)
    #
    # for ax in [ax1, ax2, ax3, ax4]:
    #     ax.set_xticks(xlabel_idx_vals)
    #     ax.set_yticks(ylabel_idx_vals)
    #     ax.set_xticklabels([r"$-2$", r"$-1$", r"$0$", r"$1$", r"$2$", ], fontsize=xticklabel_fontsize)
    #     ax.set_yticklabels([r"$0$", r"$1$", r"$2$"], fontsize=yticklabel_fontsize)
    #     ax.set_xlabel(r"$k_x$", fontsize=xlabel_fontsize)
    #     ax.set_ylabel(r"$k_y$", fontsize=ylabel_fontsize)
    #
    # for cbax in [cbax1, cbax2]:
    #     cbax.set_title(r"$\vert \tilde{\phi}_k\vert^2$", fontsize=cbax_title_fontsize)
    #
    # plt.savefig("fluxes_master.eps")
    # plt.close()

    return

def plot_properties_leapfrog_sim():
    """ """
    pickle_longname = leapfrog_pickle
    myfile = open(pickle_longname, "rb")
    [t, kx, ky, phi2, phi2_tkxky] = pickle.load(myfile)
    myfile.close()
    ## Rearrange the phi2_tkxky and kx so it goes in increasing value of kx
    sort_idxs = np.argsort(kx)
    kx = kx[sort_idxs]
    phi2_tkxky = phi2_tkxky[:,sort_idxs,:]
    print("kx = ", kx)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(t, phi2)
    ax1.scatter(t, phi2, c="red")
    for idx in [0,15,25,-1]:
        ax1.scatter(t[idx], phi2[idx], c="green")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$t$")
    ax1.set_ylabel(r"phi2")
    plt.show()

    col1_lhs = 0.1
    row1_bottom = 0.1
    plot_dim_x = 0.75
    plot_dim_y = plot_dim_x
    col1_cb_lhs = 0.9
    cb_width = 0.05

    phi2_kxky = phi2_tkxky[0,:,:]
    fig = plt.figure(figsize=[12, 12])
    ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
    cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
    cbar = make_colorplot(fig, ax1, cbax1, kx, ky, phi2_kxky, zonal=True)
    ax1.set_xlabel(r"$k_x$")
    ax1.set_ylabel(r"$k_y$")
    plt.show()

    phi2_kxky = phi2_tkxky[15,:,:]
    fig = plt.figure(figsize=[12, 12])
    ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
    cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
    cbar = make_colorplot(fig, ax1, cbax1, kx, ky, phi2_kxky, zonal=True)
    ax1.set_xlabel(r"$k_x$")
    ax1.set_ylabel(r"$k_y$")
    plt.show()

    phi2_kxky = phi2_tkxky[25,:,:]
    fig = plt.figure(figsize=[12, 12])
    ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
    cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
    cbar = make_colorplot(fig, ax1, cbax1, kx, ky, phi2_kxky, zonal=True)
    ax1.set_xlabel(r"$k_x$")
    ax1.set_ylabel(r"$k_y$")
    plt.show()

    phi2_kxky = phi2_tkxky[-1,:,:]
    fig = plt.figure(figsize=[12, 12])
    ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
    cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
    cbar = make_colorplot(fig, ax1, cbax1, kx, ky, phi2_kxky, zonal=True)
    ax1.set_xlabel(r"$k_x$")
    ax1.set_ylabel(r"$k_y$")
    plt.show()


    return

def plot_properties_nisl_sim():
    """ """
    pickle_longname = nisl_pickle
    myfile = open(pickle_longname, "rb")
    [t, kx, ky, phi2, phi2_tkxky] = pickle.load(myfile)
    myfile.close()
    ## Rearrange the phi2_tkxky and kx so it goes in increasing value of kx
    sort_idxs = np.argsort(kx)
    kx = kx[sort_idxs]
    phi2_tkxky = phi2_tkxky[:,sort_idxs,:]
    print("kx = ", kx)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(t, phi2)
    ax1.scatter(t, phi2, c="red")
    for idx in [0,15,35,-1]:
        ax1.scatter(t[idx], phi2[idx], c="green")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$t$")
    ax1.set_ylabel(r"phi2")
    plt.show()

    col1_lhs = 0.1
    row1_bottom = 0.1
    plot_dim_x = 0.75
    plot_dim_y = plot_dim_x
    col1_cb_lhs = 0.9
    cb_width = 0.05

    phi2_kxky = phi2_tkxky[0,:,:]
    fig = plt.figure(figsize=[12, 12])
    ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
    cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
    cbar = make_colorplot(fig, ax1, cbax1, kx, ky, phi2_kxky, zonal=True)
    ax1.set_xlabel(r"$k_x$")
    ax1.set_ylabel(r"$k_y$")
    plt.show()

    phi2_kxky = phi2_tkxky[15,:,:]
    fig = plt.figure(figsize=[12, 12])
    ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
    cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
    cbar = make_colorplot(fig, ax1, cbax1, kx, ky, phi2_kxky, zonal=True)
    ax1.set_xlabel(r"$k_x$")
    ax1.set_ylabel(r"$k_y$")
    plt.show()

    phi2_kxky = phi2_tkxky[35,:,:]
    fig = plt.figure(figsize=[12, 12])
    ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
    cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
    cbar = make_colorplot(fig, ax1, cbax1, kx, ky, phi2_kxky, zonal=True)
    ax1.set_xlabel(r"$k_x$")
    ax1.set_ylabel(r"$k_y$")
    plt.show()

    phi2_kxky = phi2_tkxky[-1,:,:]
    fig = plt.figure(figsize=[12, 12])
    ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
    cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
    cbar = make_colorplot(fig, ax1, cbax1, kx, ky, phi2_kxky, zonal=True)
    ax1.set_xlabel(r"$k_x$")
    ax1.set_ylabel(r"$k_y$")
    plt.show()


    return

def plot_properties_isl_sim():
    """ """
    pickle_longname = isl_pickle
    myfile = open(pickle_longname, "rb")
    [t, kx, ky, phi2, phi2_tkxky] = pickle.load(myfile)
    myfile.close()
    ## Rearrange the phi2_tkxky and kx so it goes in increasing value of kx
    sort_idxs = np.argsort(kx)
    kx = kx[sort_idxs]
    phi2_tkxky = phi2_tkxky[:,sort_idxs,:]
    print("kx = ", kx)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(t, phi2)
    ax1.scatter(t, phi2, c="red")
    for idx in [0,15,35,-1]:
        ax1.scatter(t[idx], phi2[idx], c="green")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$t$")
    ax1.set_ylabel(r"phi2")
    plt.show()

    col1_lhs = 0.1
    row1_bottom = 0.1
    plot_dim_x = 0.75
    plot_dim_y = plot_dim_x
    col1_cb_lhs = 0.9
    cb_width = 0.05

    phi2_kxky = phi2_tkxky[0,:,:]
    fig = plt.figure(figsize=[12, 12])
    ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
    cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
    cbar = make_colorplot(fig, ax1, cbax1, kx, ky, phi2_kxky, zonal=True)
    ax1.set_xlabel(r"$k_x$")
    ax1.set_ylabel(r"$k_y$")
    plt.show()

    phi2_kxky = phi2_tkxky[15,:,:]
    fig = plt.figure(figsize=[12, 12])
    ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
    cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
    cbar = make_colorplot(fig, ax1, cbax1, kx, ky, phi2_kxky, zonal=True)
    ax1.set_xlabel(r"$k_x$")
    ax1.set_ylabel(r"$k_y$")
    plt.show()

    phi2_kxky = phi2_tkxky[35,:,:]
    fig = plt.figure(figsize=[12, 12])
    ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
    cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
    cbar = make_colorplot(fig, ax1, cbax1, kx, ky, phi2_kxky, zonal=True)
    ax1.set_xlabel(r"$k_x$")
    ax1.set_ylabel(r"$k_y$")
    plt.show()

    phi2_kxky = phi2_tkxky[-1,:,:]
    fig = plt.figure(figsize=[12, 12])
    ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
    cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
    cbar = make_colorplot(fig, ax1, cbax1, kx, ky, phi2_kxky, zonal=True)
    ax1.set_xlabel(r"$k_x$")
    ax1.set_ylabel(r"$k_y$")
    plt.show()

    return

def get_data_from_pickle(pickle_longname, fluxes=False):
    """ """

    myfile = open(pickle_longname, "rb")
    if fluxes:
        [t, kx, ky, phi2, phi2_tkxky, parsed_pflx, parsed_vflx, parsed_qflx] = pickle.load(myfile)
    else:
        [t, kx, ky, phi2, phi2_tkxky] = pickle.load(myfile)

    myfile.close()

    sort_idxs = np.argsort(kx)
    kx = kx[sort_idxs]
    phi2_tkxky = phi2_tkxky[:,sort_idxs,:]

    if fluxes:
        parsed_pflx = parsed_pflx[:,:,sort_idxs,:]
        parsed_vflx = parsed_vflx[:,:,sort_idxs,:]
        parsed_qflx = parsed_qflx[:,:,sort_idxs,:]
        return [t, kx, ky, phi2, phi2_tkxky, parsed_pflx, parsed_vflx, parsed_qflx]
    else:
        return [t, kx, ky, phi2, phi2_tkxky]

def make_master_leapfrog_nisl_isl_comparison():
    """ """

    [t_isl, kx, ky, phi2_isl, phi2_tkxky_isl, parsed_pflx_isl,
     parsed_vflx_isl, parsed_qflx_isl] = get_data_from_pickle(isl_pickle, fluxes=True)
    [t_nisl, kx, ky, phi2_nisl, phi2_tkxky_nisl, parsed_pflx_nisl,
     parsed_vflx_nisl, parsed_qflx_nisl] = get_data_from_pickle(nisl_pickle, fluxes=True)
    [t_leapfrog, kx, ky, phi2_leapfrog, phi2_tkxky_leapfrog, parsed_pflx_leapfrog,
     parsed_vflx_leapfrog, parsed_qflx_leapfrog] = get_data_from_pickle(leapfrog_pickle, fluxes=True)
    [t_master, kx, ky, phi2_master, phi2_tkxky_master, parsed_pflx_master,
     parsed_vflx_master, parsed_qflx_master] = get_data_from_pickle(master_cfl025_pickle, fluxes=True)

    marker_size = 5
    linewidth=2
    legend_fontsize = 16
    xlabel_fontsize = 30
    ylabel_fontsize = 30
    xticklabel_fontsize = 18
    yticklabel_fontsize = 18

    left = 0.16
    right =0.95
    top = 0.95
    bottom = 0.14
    width = right - left
    height = top - bottom
    fig = plt.figure(figsize=(12,10))
    ax1 = fig.add_axes((left, bottom, width, height))
    linewidths = [5, 5, 4, 3, 2, 3]
    linestyles = ["-", "--", "-.", (0,(3,1,2,1)), "-", "--", "-."]
    col_counters = [0, 3, 4, 6]
    t_list = [t_master, t_leapfrog, t_nisl, t_isl]
    phi2_list = [phi2_master, phi2_leapfrog, phi2_nisl, phi2_isl]
    labels=[r"SSP RK2", r"Leapfrog", r"NISL", r"SL (bilinear interpolation)"]

    for counter in [0, 1, 2, 3]:
        col_counter = col_counters[counter]
        ax1.plot(t_list[counter], phi2_list[counter], lw=linewidths[counter],
                 label=labels[counter], ls=linestyles[col_counter], c=default_cmap(col_counter))

    ax1.grid(True)
    ax1.legend(loc="lower right", fontsize=legend_fontsize)
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$\tilde{t}$", fontsize=xlabel_fontsize)
    ax1.set_ylabel(r"$ \sum_{\tilde{k}_x, \tilde{k}_y} \vert \tilde{\phi}_k\vert^2$", fontsize=ylabel_fontsize)
    ax1.set_xlim(0, 400)
    ax1.set_ylim(2E-4, 1E6)
    ax1.set_xticks([0, 100, 200, 300, 400])
    ax1.set_xticklabels([r"$0$", r"$100$", r"$200$", r"$300$", r"$400$"], fontsize=xticklabel_fontsize)
    ax1.set_yticks([0.01, 1, 100, 1E4, 1E6])
    ax1.set_yticklabels([r"$10^{-2}$", r"$10^0$", r"$10^2$", r"$10^4$", r"$10^6$"], fontsize=yticklabel_fontsize)
    plt.savefig("nonlinear_sim_comparisons_phi2_t.eps")
    plt.close()



    ########### A look at phi2(kx, ky) ##################


    cbarticklabel_fontsize = 20
    title_fontsize = 36
    ## imshow makes the plot, but in index-space rather than (kx, ky)
    # To have meaningful labels, we want to find indexes corresponding to the
    # desired values of kx, ky
    kx_label_vals = np.array([-2, -1, 0, 1, 2])
    ky_label_vals = np.array([0, 1, 2])
    # idx space goes [0, 1, 2, . . . , nkx(y)-1]
    # kx(y) space goes [kx(y)_min,  . . . , kx(y)_max]
    # so kx = ((kx_max - kx_min)/(nkx-1))*kxidx + kx_min
    # so kxidx = (kx - kx_min) * ((nkx-1)/(kx_max - kx_min))
    kx_min = np.min(kx) ; kx_max = np.max(kx)
    ky_min = np.min(ky) ; ky_max = np.max(ky)

    kx_extent = kx_max - kx_min
    ky_extent = ky_max - ky_min
    nkx = len(kx) ; nky = len(ky)
    xlabel_idx_vals = (kx_label_vals - kx_min) * ((nkx-1)/kx_extent)
    ylabel_idx_vals = (ky_label_vals - ky_min) * ((nky-1)/ky_extent)

    fig = plt.figure(figsize=[12, 12])

    left = 0.1
    right = 0.92
    top = 0.94
    bottom = 0.08

    # row1_bottom = 0.1
    subplot_hspace = 0.1
    subplot_vspace = 0.1
    subplot_width = (right - left - subplot_hspace)/2
    subplot_height = (top - bottom - subplot_vspace)/2
    cb_width = 0.02
    cb_hspace = 0.01
    plot_width = subplot_width - cb_width - cb_hspace
    plot_height = subplot_height
    print("plot_width, plot_height = ", plot_width, plot_height)

    ax_col1_lhs = left
    ax_col2_lhs = left + subplot_width + subplot_hspace
    cbax_col1_lhs = left + plot_width + cb_hspace
    cbax_col2_lhs = left + subplot_width + subplot_hspace + plot_width + cb_hspace
    row1_bottom = bottom + subplot_vspace + subplot_height

    xticklabel_fontsize = 24
    yticklabel_fontsize = 24
    xlabel_fontsize = 30
    ylabel_fontsize = 30
    cbax_title_fontsize = 30
    # plot_dim_x = 0.75
    # plot_dim_y = plot_dim_x
    # col1_cb_lhs = 0.9

    ax1 = fig.add_axes([ax_col1_lhs, row1_bottom, plot_width, plot_height])
    cbax1 = fig.add_axes([cbax_col1_lhs, row1_bottom, cb_width, plot_height])
    ax2 = fig.add_axes([ax_col2_lhs, row1_bottom, plot_width, plot_height])
    cbax2 = fig.add_axes([cbax_col2_lhs, row1_bottom, cb_width, plot_height])
    ax3 = fig.add_axes([ax_col1_lhs, bottom, plot_width, plot_height])
    cbax3 = fig.add_axes([cbax_col1_lhs, bottom, cb_width, plot_height])
    ax4 = fig.add_axes([ax_col2_lhs, bottom, plot_width, plot_height])
    cbax4 = fig.add_axes([cbax_col2_lhs, bottom, cb_width, plot_height])

    axes = [ax1, ax2, ax3, ax4]
    cbaxes = [cbax1, cbax2, cbax3, cbax4]

    cbar1 = make_colorplot(fig, ax1, cbax1, kx, ky, np.log10(phi2_tkxky_master[-1,:,:]), zonal=True, log=False)
    cbar2 = make_colorplot(fig, ax2, cbax2, kx, ky, np.log10(phi2_tkxky_leapfrog[-1,:,:]), zonal=True, log=False)
    cbar3 = make_colorplot(fig, ax3, cbax3, kx, ky, np.log10(phi2_tkxky_nisl[-1,:,:]), zonal=True, log=False)
    cbar4 = make_colorplot(fig, ax4, cbax4, kx, ky, np.log10(phi2_tkxky_isl[-1,:,:]), zonal=True, log=False)

    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_xticks(xlabel_idx_vals)
        ax.set_yticks(ylabel_idx_vals)
        ax.set_xticklabels([r"$-2$", r"$-1$", r"$0$", r"$1$", r"$2$", ], fontsize=xticklabel_fontsize)
        ax.set_yticklabels([r"$0$", r"$1$", r"$2$"], fontsize=yticklabel_fontsize)

    for ax in [ax1, ax3]:
        ax.set_ylabel(r"$\tilde{k}_y$", fontsize=ylabel_fontsize)

    for ax in [ax3, ax4]:
        ax.set_xlabel(r"$\tilde{k}_x$", fontsize=xlabel_fontsize)

    for cbax in [cbax1, cbax2]:
        cbax.set_title(r"$\vert \tilde{\phi}_k\vert^2$", fontsize=cbax_title_fontsize)

    cbar1.set_ticks([-3, -2, -1, 0])
    cbax1.set_yticklabels([r"$10^{-3}$", r"$10^{-2}$", r"$10^{-1}$", r"$10^0$"], fontsize=cbarticklabel_fontsize)
    cbar2.set_ticks([1, 2, 3])
    cbax2.set_yticklabels([r"$10^1$", r"$10^2$", r"$10^3$"], fontsize=cbarticklabel_fontsize)
    cbar3.set_ticks([-4, -2, 0, 2, 4])
    cbax3.set_yticklabels([r"$10^{-4}$", r"$10^{-2}$", r"$10^{0}$", r"$10^2$", r"$10^4$"], fontsize=cbarticklabel_fontsize)
    cbar4.set_ticks([-4, 0, 4])
    cbax4.set_yticklabels([r"$10^{-4}$", r"$10^{0}$", r"$10^4$"], fontsize=cbarticklabel_fontsize)

    ax1.set_title(r"RK2", fontsize=title_fontsize)
    ax2.set_title(r"Leapfrog", fontsize=title_fontsize)
    ax3.set_title(r"NISL", fontsize=title_fontsize)
    ax4.set_title(r"SL (bilin. interp.)", fontsize=title_fontsize)
    plt.savefig("phi2_different_sims.eps")
    plt.close()


    ########### A look at heat fluxes ##################

    ## imshow makes the plot, but in index-space rather than (kx, ky)
    # To have meaningful labels, we want to find indexes corresponding to the
    # desired values of kx, ky
    kx_label_vals = np.array([-2, -1, 0, 1, 2])
    ky_label_vals = np.array([0, 1, 2])
    # idx space goes [0, 1, 2, . . . , nkx(y)-1]
    # kx(y) space goes [kx(y)_min,  . . . , kx(y)_max]
    # so kx = ((kx_max - kx_min)/(nkx-1))*kxidx + kx_min
    # so kxidx = (kx - kx_min) * ((nkx-1)/(kx_max - kx_min))
    kx_min = np.min(kx) ; kx_max = np.max(kx)
    ky_min = np.min(ky) ; ky_max = np.max(ky)

    kx_extent = kx_max - kx_min
    ky_extent = ky_max - ky_min
    nkx = len(kx) ; nky = len(ky)
    xlabel_idx_vals = (kx_label_vals - kx_min) * ((nkx-1)/kx_extent)
    ylabel_idx_vals = (ky_label_vals - ky_min) * ((nky-1)/ky_extent)

    fig = plt.figure(figsize=[12, 12])

    left = 0.1
    right = 0.95
    top = 0.92
    bottom = 0.1

    # row1_bottom = 0.1
    subplot_hspace = 0.1
    subplot_vspace = 0.1
    subplot_width = (right - left - subplot_hspace)/2
    subplot_height = (top - bottom - subplot_vspace)/2
    cb_width = 0.02
    cb_hspace = 0.02
    plot_width = subplot_width - cb_width - cb_hspace
    plot_height = subplot_height
    print("plot_width, plot_height = ", plot_width, plot_height)

    ax_col1_lhs = left
    ax_col2_lhs = left + subplot_width + subplot_hspace
    cbax_col1_lhs = left + plot_width + cb_hspace
    cbax_col2_lhs = left + subplot_width + subplot_hspace + plot_width + cb_hspace
    row1_bottom = bottom + subplot_vspace + subplot_height

    xticklabel_fontsize = 24
    yticklabel_fontsize = 24
    xlabel_fontsize = 30
    ylabel_fontsize = 30
    cbax_title_fontsize = 30
    # plot_dim_x = 0.75
    # plot_dim_y = plot_dim_x
    # col1_cb_lhs = 0.9

    ax1 = fig.add_axes([ax_col1_lhs, row1_bottom, plot_width, plot_height])
    cbax1 = fig.add_axes([cbax_col1_lhs, row1_bottom, cb_width, plot_height])
    ax2 = fig.add_axes([ax_col2_lhs, row1_bottom, plot_width, plot_height])
    cbax2 = fig.add_axes([cbax_col2_lhs, row1_bottom, cb_width, plot_height])
    ax3 = fig.add_axes([ax_col1_lhs, bottom, plot_width, plot_height])
    cbax3 = fig.add_axes([cbax_col1_lhs, bottom, cb_width, plot_height])
    ax4 = fig.add_axes([ax_col2_lhs, bottom, plot_width, plot_height])
    cbax4 = fig.add_axes([cbax_col2_lhs, bottom, cb_width, plot_height])

    axes = [ax1, ax2, ax3, ax4]
    cbaxes = [cbax1, cbax2, cbax3, cbax4]

    cbar = make_colorplot(fig, ax1, cbax1, kx, ky, parsed_qflx_master[-1,0,:,:], zonal=True, log=False)
    cbar = make_colorplot(fig, ax2, cbax2, kx, ky, parsed_qflx_leapfrog[-1,0,:,:], zonal=True, log=False)
    cbar = make_colorplot(fig, ax3, cbax3, kx, ky, parsed_qflx_nisl[-1,0,:,:], zonal=True, log=False)
    cbar = make_colorplot(fig, ax4, cbax4, kx, ky, parsed_qflx_isl[-1,0,:,:], zonal=True, log=False)

    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_xticks(xlabel_idx_vals)
        ax.set_yticks(ylabel_idx_vals)
        ax.set_xticklabels([r"$-2$", r"$-1$", r"$0$", r"$1$", r"$2$", ], fontsize=xticklabel_fontsize)
        ax.set_yticklabels([r"$0$", r"$1$", r"$2$"], fontsize=yticklabel_fontsize)
        ax.set_xlabel(r"$k_x$", fontsize=xlabel_fontsize)
        ax.set_ylabel(r"$k_y$", fontsize=ylabel_fontsize)

    for cbax in [cbax1, cbax2]:
        cbax.set_title(r"Heat flux", fontsize=cbax_title_fontsize)

    plt.savefig("heat_flux_different_sims.eps")
    plt.close()

    return

if __name__ == "__main__":
    print("Hello world")
    # plot_properties_master_sim(cfl_cushion="0.25")
    # plot_properties_leapfrog_sim()
    # plot_properties_nisl_sim()
    # plot_properties_isl_sim()
    # make_master_cushion_comparison()
    make_master_leapfrog_nisl_isl_comparison()
