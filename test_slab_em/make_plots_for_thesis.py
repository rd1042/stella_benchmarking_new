"""Make plots for electromagnetic slab sims"""

import sys
sys.path.append("../postprocessing_tools")
from plotting_helper import make_beta_scan_plots, make_comparison_plots, plot_gmvus, plot_gzvs
from helper_ncdf_new import view_ncdf_variables, extract_data_from_ncdf, extract_data_from_ncdf_with_xarray, view_ncdf_variables_with_xarray
from extract_sim_data import find_bpar_phi_ratio, find_apar_phi_ratio
import matplotlib.pyplot as plt
import numpy as np
import glob
import re
import pickle
# import ast

IMAGE_DIR = "./images/"

# basic_em_sim = "stella_sims/input_slab_ky0.1_explicit"
# mandell_beta1_kperp1 = "stella_sims/input_slab_ky0.1_explicit_mandell1"
# mandell_beta1_kperp001_new = "mandell_sims/input_slab_ky0.01_beta1_new"
# mandell_beta1_kperp1_new = "mandell_sims/input_slab_ky1_beta1_new"
# mandell_sf1_kperp1_new = "mandell_sims/input_slab_ky1_sf1_new"
# mandell_sf1_kperp001_new = "mandell_sims/input_slab_ky0.01_sf1_new"
# mandell_sf0002_kperp001_new = "mandell_sims/input_slab_ky0.01_sf0.002_new"
# mandell_sf00002_kperp001_new = "mandell_sims/input_slab_ky0.01_sf0.0002_new"
# mandell_sf000002_kperp001_new = "mandell_sims/input_slab_ky0.01_sf0.00002_new"
# mandell_sf10_kperp001_new = "mandell_sims/input_slab_ky0.01_sf10_new"
# mandell_beta1_kperp1_long_t = "stella_sims/input_slab_ky0.1_explicit_mandell2"
# mandell_beta1_kperp1_long_t_marconi = "mandell_sims/input_slab_ky0.1_beta1"


def compare_sims_for_thesis(field_name, outnc_longnames, sim_types, labels, save_name="a.eps",
            ax1_yticks=None, ax1_yticklabels=None, ax2_yticks=None, ax2_yticklabels=None):
    """ """

    ylabel_fontsize = 30
    xlabel_fontsize = 30
    legend_fontsize = 18
    linestyles=((0,(1,0)), (0, (2,2)), (0,(5,1,1,1)), (0,(1,1)),
                (0, (3,2,2,3)), (0, (1,0)) )
    # linewidths = [6, 4, 3, 3, 3, 2]
    # linestyles=((0,(1,0)), (0, (2,2)), (0,(5,1,1,1)), (0,(1,1)), (0,(1,0)), (0, (2,2)))
    # linewidths = [4, 4, 3, 3, 2, 2]
    y_ticklabelfontsize = 25
    x_ticklabelfontsize = y_ticklabelfontsize
    my_xticklength = 5
    my_xtickwidth = 1
    my_yticklength = 5
    my_ytickwidth = 1

    default_cols = plt.get_cmap("tab10")

    left = 0.15
    right = 0.985
    top = 0.985
    bottom = 0.1
    vspace=0.1
    height = (top - bottom - vspace)/2
    width = right - left
    my_linewidth = 3

    fig = plt.figure(figsize=(12,12))
    ax1 = fig.add_axes((left, bottom + height + vspace, width, height))
    ax2 = fig.add_axes((left, bottom, width, height))

    t_list = []
    z_list = []
    phi_vs_t_list = []

    for idx, outnc_longname in enumerate(outnc_longnames):
        if sim_types[idx] == "stella":
            t, z, phi_vs_t, beta = extract_data_from_ncdf_with_xarray(outnc_longname, "t", 'zed', 'phi_vs_t', 'beta')
            phi_vs_t = np.array(phi_vs_t[:,0,:,0,0,:])
            # shape is now (time, zed, ri)
        elif sim_types[idx] == "gs2":
            t, z, phi_vs_t, beta = extract_data_from_ncdf_with_xarray(outnc_longname, "t", 'theta', 'phi_t', 'beta')
            phi_vs_t = np.array(phi_vs_t[:,0,0,:,:])
            # shape is now (time, zed, ri)
        else:
            print("sim type not recognised! Aborting")
            sys.exit()
        t_list.append(t); z_list.append(z); phi_vs_t_list.append(phi_vs_t);
    for idx, outnc_longname in enumerate(outnc_longnames):
        t = t_list[idx]; z = z_list[idx]; phi_vs_t = phi_vs_t_list[idx];
        zval= float(z[int(len(z)*0.5)])
        if abs(zval) > 1E-5:
            print("Error! zval!=0 . zval = ", zval)
            sys.exit()
        phi_t = phi_vs_t[:,int(len(z)*0.5),0]
        if idx == 0:
            phi_t_orig = phi_t
            z_orig = zval
            t_orig = t
            phi_t_normalising_factor = np.max(abs(phi_t_orig))
        ax1.plot(t, phi_t, label=(labels[idx]), ls=linestyles[idx], lw=my_linewidth)


        if len(t_orig) == len(t):
            if np.max(abs(t - t_orig)) < 1E-5:
                if idx != 0:    # Don't plot the first, as this is what we're comparing to
                    ax2.plot(t, (abs((phi_t - phi_t_orig)/phi_t_normalising_factor)*100),
                                        c=default_cols(idx), ls=linestyles[idx], lw=my_linewidth)
            else:
                print("t vals don't match")
                print("t_orig - tval = ", t_orig - t)
        else:
            print("t lengths not equal")
            print("outnc_longname = ", outnc_longname)

    if field_name == "phi":
        field_ylabel_ax1 = r"$\tilde{\phi}_k$"
        field_ylabel_ax2 = r"$\tilde{\phi}_k$ error (%)"
    if field_name == "apar":
        field_ylabel_ax1 = r"$\tilde{A}_{\parallel, k}$"
        field_ylabel_ax2 = r"$\tilde{A}_{\parallel, k}$ error (%)"
    if field_name == "bpar":
        field_ylabel_ax1 = r"$\tilde{B}_{\parallel, k}$"
        field_ylabel_ax2 = r"$\tilde{B}_{\parallel, k}$ error (%)"

    # ax1.set_xlim((-1, 1))
    # ax2.set_xlim((-1, 1))
    # ax1.set_xlim((0,20))
    # ax1.set_ylim((-6, 6))
    if ax1_yticks is not None:
        ax1.set_yticks(ax1_yticks)
        ax1.set_yticklabels(ax1_yticklabels, fontsize=y_ticklabelfontsize)
    if ax2_yticks is not None:
        ax2.set_yticks(ax2_yticks)
        ax2.set_yticklabels(ax2_yticklabels, fontsize=y_ticklabelfontsize)
    ax1.tick_params("x", length=my_xticklength, width=my_xtickwidth, direction="out")
    ax2.tick_params("x", length=my_xticklength, width=my_xtickwidth, direction="out")
    ax1.tick_params("y", length=my_yticklength, width=my_ytickwidth, direction="out")
    # ax1.set_xticks([-1, 0, 1])
    # ax1.set_xticklabels([])
    # ax2.set_xticks([-1, 0, 1])
    # ax2.set_xticklabels(["-1", "0", "1"], fontsize=x_ticklabelfontsize)
    ax1.legend(loc="best", fontsize=legend_fontsize)
    ax1.set_ylabel(field_ylabel_ax1, fontsize=ylabel_fontsize)
    ax2.set_ylabel(field_ylabel_ax2, fontsize=ylabel_fontsize, labelpad=20)
    ax2.set_xlabel(r"$\tilde{t}$", fontsize=xlabel_fontsize)
    ax1.set_xlabel(r"$\tilde{t}$", fontsize=xlabel_fontsize)

    # plt.show()
    plt.savefig(save_name)
    plt.close()


    return

# def compare_sims_with_apar(outnc_longnames, labels):
#     """ """
#
#     fig = plt.figure()
#     ax1 = fig.add_subplot(221)
#     ax2 = fig.add_subplot(223, sharex=ax1)
#     ax3 = fig.add_subplot(222)
#     ax4 = fig.add_subplot(224, sharex=ax3)
#     for idx, outnc_longname in enumerate(outnc_longnames):
#         t, z, phi_vs_t, apar_vs_t, beta = extract_data_from_ncdf(outnc_longname, "t", 'zed', 'phi_vs_t', 'apar_vs_t', 'beta')
#         phi_vs_t = phi_vs_t[:,0,:,0,0]
#         apar_vs_t = apar_vs_t[:,0,:,0,0]
#         print("len(t) = ", len(t))
#         print("len(z) = ", len(z))
#         print("phi_vs_t.shape = ", phi_vs_t.shape)
#         # print("phi_vs_t = ", phi_vs_t)
#         # print("phi_vs_t.imag = ", phi_vs_t.imag)
#         # sys.exit()
#         ax1.plot(z, phi_vs_t[0,:].real, label=(labels[idx] + ", t=0"))
#         ax2.plot(z, phi_vs_t[0,:].imag)
#         ax3.plot(z, apar_vs_t[0,:].real, label=(labels[idx] + ", t=0" ))
#         ax4.plot(z, apar_vs_t[0,:].imag)
#     ax2.set_xlabel("z")
#     ax4.set_xlabel("z")
#     ax2.set_ylabel("Im(phi)")
#     ax1.set_ylabel("Re(phi)")
#     ax4.set_ylabel("Im(apar)")
#     ax3.set_ylabel("Re(apar)")
#     ax1.legend(loc="best")
#     ax3.legend(loc="best")
#     for ax in [ax1, ax2, ax3, ax4]:
#         ax.grid(True)
#     plt.show()
#
#     fig = plt.figure()
#     ax1 = fig.add_subplot(221)
#     ax2 = fig.add_subplot(223, sharex=ax1)
#     ax3 = fig.add_subplot(222)
#     ax4 = fig.add_subplot(224, sharex=ax3)
#     for idx, outnc_longname in enumerate(outnc_longnames):
#         t, z, phi_vs_t, apar_vs_t, beta = extract_data_from_ncdf(outnc_longname, "t", 'zed', 'phi_vs_t', 'apar_vs_t', 'beta')
#         phi_vs_t = phi_vs_t[:,0,:,0,0]
#         apar_vs_t = apar_vs_t[:,0,:,0,0]
#         print("len(t) = ", len(t))
#         print("len(z) = ", len(z))
#         print("phi_vs_t.shape = ", phi_vs_t.shape)
#         # print("phi_vs_t = ", phi_vs_t)
#         # print("phi_vs_t.imag = ", phi_vs_t.imag)
#         # sys.exit()
#         ax1.plot(t, phi_vs_t[:,int(len(z)*0.5)].real, label=(labels[idx] + ", z=" + str(z[int(len(z)*0.5)]) ))
#         ax2.plot(t, phi_vs_t[:,int(len(z)*0.5)].imag)
#         ax3.plot(t, apar_vs_t[:,int(len(z)*0.25)].real, label=(labels[idx] + ", z=" + str(z[int(len(z)*0.25)]) ))
#         ax4.plot(t, apar_vs_t[:,int(len(z)*0.25)].imag)
#     ax2.set_xlabel("t")
#     ax4.set_xlabel("t")
#     ax2.set_ylabel("Im(phi)")
#     ax1.set_ylabel("Re(phi)")
#     ax4.set_ylabel("Im(apar)")
#     ax3.set_ylabel("Re(apar)")
#     ax1.legend(loc="best")
#     ax3.legend(loc="best")
#     for ax in [ax1, ax2, ax3, ax4]:
#         ax.grid(True)
#     plt.show()
#
#     return

def benchmark_stella_vs_gs2_fapar0_fbpar0_for_thesis():
    """ """

    gs2_implicit_dt1Em2 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar0_fbpar0/input_nz36_dt1E-2.out.nc"
    gs2_implicit_dt5Em4 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar0_fbpar0/input_nz36_dt5E-4.out.nc"
    gs2_implicit_dt25Em5 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar0_fbpar0/input_nz36_dt25E-5.out.nc"
    stella_implicit_dt1Em2 = "sims/for_thesis_stella_ky1_beta1_zero_gradients_fapar0_fbpar0/input_implicit_nz36_dt1E-2.out.nc"
    stella_implicit_dt5Em4 = "sims/for_thesis_stella_ky1_beta1_zero_gradients_fapar0_fbpar0/input_implicit_nz36_dt5E-4.out.nc"
    # stella_implicit_dt25Em5 = "sims/for_thesis_stella_ky1_beta1_zero_gradients_fapar0_fbpar0/input_implicit_nz36_dt25E-5.out.nc"

    compare_sims_for_thesis("phi",
            [gs2_implicit_dt25Em5, gs2_implicit_dt5Em4, gs2_implicit_dt1Em2,
             stella_implicit_dt1Em2, stella_implicit_dt5Em4, # stella_implicit_dt25Em5
             ],
            ["gs2", "gs2", "gs2",
             "stella", "stella", # "stella"
             ],
            [
            "gs2(dt=2.5E-5)", "gs2(dt=5E-4)", "gs2(dt=1E-2)",
            "stella (dt=1E-2)", "stella (dt=5E-4)",  # "stella (dt=5E-2)",
            ],
             save_name="for_thesis_phi_t_fapar0_fbpar0_dt_scan.eps")

    return

def benchmark_stella_vs_gs2_fapar1_fbpar0_for_thesis():
    """ """
    # gs2_sim = "sims/gs2_ky1_beta1_zero_gradients_fapar1_fbpar0/input_long.out.nc"
    # gs2_sim_tup0 = "sims/gs2_ky1_beta1_zero_gradients_fapar1_fbpar0/input_long_fexpr_0.5.out.nc"
    # gs2_sim_tup0_zup0 = "sims/gs2_ky1_beta1_zero_gradients_fapar1_fbpar0/input_long_fexpr_0.5_bakdif0.out.nc"
    # stella_sim = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_fully_explicit.out.nc"
    # stella_sim_str_impl = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_implicit_streaming.out.nc"
    # stella_sim_str_impl_zup0_tup0 = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_implicit_streaming_zup0_tup0.out.nc"
    # stella_sim_str_impl_zup0 = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_implicit_streaming_zup0.out.nc"
    # stella_sim_str_impl_tup0 = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_implicit_streaming_tup0.out.nc"
    # compare_sims_for_thesis("phi", [gs2_sim, stella_sim_str_impl,
    #               stella_sim_str_impl_zup0, stella_sim_str_impl_tup0,
    #               ],
    #              ["gs2",   "stella",
    #              "stella", "stella",],
    #              ["gs2", "stella (implicit, zupw=tupw=0.02)",
    #               "stella (implicit, zupw=0, tupw=0.02)", "stella (implicit, zupw=0.2, tupw=0)",
    #               ],
    #              save_name="phi_t_fapar1_fbpar0.eps")
    #
    #
    # stella_sim_new_init = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar1/input_fully_explicit_init_matching_gs2.out.nc"
    # stella_sim_new_init_ntheta144 = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar1/input_fully_explicit_init_matching_gs2_ntheta144.out.nc"
    # gs2_long_sim = "sims/gs2_ky1_beta1_zero_gradients_fapar1_fbpar1/input_long.out.nc"
    # compare_sims_for_thesis("phi", [gs2_long_sim, stella_sim_new_init, stella_sim_new_init_ntheta144],
    #         ["gs2", "stella", "stella"],
    #         ["gs2", "stella (explicit)", "stella (explicit) nzed=144"],
    #          save_name="phi_t_fapar1_fbpar1_explicit.eps")

    gs2_implicit_dt5Em2 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar0/input_nz36_dt5E-2.out.nc"
    gs2_implicit_dt1Em2 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar0/input_nz36_dt1E-2.out.nc"
    gs2_implicit_dt5Em4 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar0/input_nz36_dt5E-4.out.nc"
    #gs2_implicit_dt25Em5 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar0/input_nz36_dt25E-5.out.nc"
    stella_implicit_dt5Em2 = "sims/for_thesis_stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_implicit_nz36_dt5E-2.out.nc"
    stella_implicit_dt1Em2 = "sims/for_thesis_stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_implicit_nz36_dt1E-2.out.nc"
    stella_implicit_dt5Em4 = "sims/for_thesis_stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_implicit_nz36_dt5E-4.out.nc"
    #stella_implicit_dt25Em5 = "sims/for_thesis_stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_implicit_nz36_dt25E-5.out.nc"

    ## explicit sims
    stella_explicit = "sims/for_thesis_stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_explicit_nz36_dt25E-5.out.nc"
    stella_explicit_centered_dgdz = "sims/for_thesis_stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_explicit_nz36_dt25E-5_centered_dgdz.out.nc"

    compare_sims_for_thesis("phi",
            [# gs2_implicit_dt25Em5,
             gs2_implicit_dt5Em4, gs2_implicit_dt1Em2, gs2_implicit_dt5Em2,
             stella_implicit_dt5Em2, stella_implicit_dt1Em2, stella_implicit_dt5Em4,
             #stella_implicit_dt25Em5
             ],
            ["gs2", "gs2", "gs2", #"gs2",
             "stella", "stella", "stella", # "stella"
             ],
            [
            # "gs2(dt=2.5E-4)",
            "gs2(dt=5E-4)", "gs2(dt=1E-2)", "gs2(dt=5E-2)",
            "stella (dt=5E-2)", "stella (dt=1E-2)", "stella (dt=5E-4)",
            # "stella (dt=2.5E-4)",
            ],
             save_name="for_thesis_phi_t_fapar1_fbpar0_dt_scan.eps")

    compare_sims_for_thesis("phi",
            [# gs2_implicit_dt25Em5,
             gs2_implicit_dt5Em4, stella_implicit_dt5Em4, stella_explicit,
             stella_explicit_centered_dgdz,
             #stella_implicit_dt25Em5
             ],
            ["gs2",
             "stella", "stella", "stella", # "stella"
             ],
            [
            # "gs2(dt=2.5E-4)",
            "gs2(dt=5E-4)",
            "stella (implicit)", "stella (explicit)", "stella (explicit, centered dg/dz)",
            # "stella (dt=2.5E-4)",
            ],
             save_name="for_thesis_phi_t_fapar1_fbpar0_explicit_centered_dgdz.eps")

    return

def benchmark_stella_src_h_stella_vs_gs2_fapar1_fbpar1():
    """ """
    gs2_implicit_dt5Em2 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar1/input_nz36_dt5E-2.out.nc"
    # stella_implicit_dt5Em2 = "sims/for_thesis_stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_implicit_nz36_dt5E-2.out.nc"
    stella_explicit_src_h_dt4Em4 = "sims/stella_comparing_implex_gbar_h_src/input_explicit_src_h.out.nc"
    stella_implicit_src_h_dt4Em4 = "sims/stella_comparing_implex_gbar_h_src/input_implicit_src_h.out.nc"
    stella_implicit_src_h_dt4Em2 = "sims/stella_comparing_implex_gbar_h_src/input_implicit_src_h_dt4E-2.out.nc"
    stella_implicit_src_h_dt4Em2_tupw002 = "sims/stella_comparing_implex_gbar_h_src/input_implicit_src_h_dt4E-2_tupw002.out.nc"

    compare_sims_for_thesis("phi",
            [# gs2_implicit_dt25Em5,
             gs2_implicit_dt5Em2,
             stella_explicit_src_h_dt4Em4,
             stella_implicit_src_h_dt4Em4,
             stella_implicit_src_h_dt4Em2,
             #stella_implicit_src_h_dt4Em2_tupw002
             #stella_implicit_dt25Em5
             ],
            ["gs2",
             "stella", "stella", "stella", #"stella"
             ],
            [
            # "gs2(dt=2.5E-4)",
            "gs2(dt=5E-4)",
            "stella (explicit, src. h)",
            "stella (impl, src. h, dt=4E-4, zupw=tupw=0)",
            "stella (impl, src. h, dt=4E-2, zupw=tupw=0)",#, "stella (explicit, centered dg/dz)",
            #"stella (impl, src. h, dt=4E-2, tupw=0.02)"#, "stella (explicit, centered dg/dz)",
            # "stella (dt=2.5E-4)",
            ],
             save_name="for_thesis_phi_t_fapar1_fbpar1_src_h.eps")

def benchamrk_stella_src_h_stella_vs_gs2_fapar1_fbpar0():
    """ """
    gs2_implicit_dt5Em2 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar0/input_nz36_dt5E-2.out.nc"
    stella_implicit_dt5Em2 = "sims/for_thesis_stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_implicit_nz36_dt5E-2.out.nc"
    stella_explicit_src_h_dt4Em4 = "sims/stella_comparing_implex_gbar_h_src/input_explicit_src_h_fapar1_fbpar0.out.nc"
    # stella_implicit_src_h_dt4Em4 = "sims/stella_comparing_implex_gbar_h_src/input_implicit_src_h.out.nc"
    stella_implicit_src_h_dt4Em2 = "sims/stella_comparing_implex_gbar_h_src/input_implicit_src_h_fapar1_fbpar0_dt4E-2.out.nc"


    compare_sims_for_thesis("phi",
            [# gs2_implicit_dt25Em5,
             gs2_implicit_dt5Em2, stella_implicit_dt5Em2,
             stella_explicit_src_h_dt4Em4,
             #stella_implicit_src_h_dt4Em4,
             stella_implicit_src_h_dt4Em2
             #stella_implicit_dt25Em5
             ],
            ["gs2",
             "stella", "stella", "stella", # "stella"
             ],
            [
            # "gs2(dt=2.5E-4)",
            "gs2(dt=5E-4)",
            "stella (implicit, dt=4E-2)",
            "stella (explicit, src. h)",
            #"stella (impl, src. h, dt=4E-4)",
            "stella (impl, src. h, dt=4E-2)"#, "stella (explicit, centered dg/dz)",
            # "stella (dt=2.5E-4)",
            ],
             save_name="for_thesis_phi_t_fapar1_fbpar0_src_h.eps")

    return

if __name__ == "__main__":
    print("Hello world")

    # benchmark_stella_vs_gs2_fapar0_fbpar0_for_thesis()
    # benchmark_stella_vs_gs2_fapar1_fbpar0_for_thesis()
    benchmark_stella_src_h_stella_vs_gs2_fapar1_fbpar1()
    # benchamrk_stella_src_h_stella_vs_gs2_fapar1_fbpar0()