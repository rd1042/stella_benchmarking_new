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
from scipy.optimize import curve_fit

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


def damped_oscillation(tdata, amp0, phase, freq, gamma, offset):
    return (offset + amp0*np.sin(freq*tdata+phase)*np.exp(gamma*tdata))

def fit_ksaw(t, phi_t, make_plot=False):
    """Fit a straight line to the natural logarithm of phi**2. If we assume that
    phi**2 is described by phi**2 = A*exp(Bt), then log(phi**2) = ln(A) + Bt"""

    guess_amp = 4.1
    guess_phase = 1.5*np.pi
    guess_freq = 1.1
    guess_gamma = -0.03
    guess_offset = 0
    initial_guesses = np.array([guess_amp,
                                guess_phase,
                                guess_freq,
                                guess_gamma,
                                guess_offset])
    # if make_plot:
    #     tmin = np.min(t) ; tmax = np.max(t)
    #     upsampled_t = np.linspace(tmin, tmax, 1000)
    #     fig = plt.figure()
    #     ax1 = fig.add_subplot(111)
    #
    #     ax1.scatter(t, phi_t, c="black", marker="x", s=10)
    #     plt.show()


    popt, pcov = curve_fit(damped_oscillation, t, phi_t, p0=initial_guesses)
    [amp0_opt, phase_opt, freq_opt, gamma_opt, offset_opt] = popt
    amp0_err = (np.sqrt(pcov[0, 0]))
    phase_err = (np.sqrt(pcov[1, 1]))
    freq_err = (np.sqrt(pcov[2, 2]))
    gamma_err = (np.sqrt(pcov[3, 3]))
    offset_err = (np.sqrt(pcov[4, 4]))

    if make_plot:
        tmin = np.min(t) ; tmax = np.max(t)
        upsampled_t = np.linspace(tmin, tmax, 1000)
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(upsampled_t, damped_oscillation(upsampled_t, guess_amp, guess_phase, guess_freq, guess_gamma, guess_offset), c="blue")
        ax1.plot(upsampled_t, damped_oscillation(upsampled_t, amp0_opt, phase_opt, freq_opt, gamma_opt, offset_opt), c="red")
        ax1.scatter(t, phi_t, c="black", marker="x", s=10)
        plt.show()
    return ([amp0_opt, phase_opt, freq_opt, gamma_opt, offset_opt],
           [amp0_err, phase_err, freq_err, gamma_err, offset_err])

def fit_frequency_and_damping_rate(outnc_longname, sim_type, **kwargs):
    """ """
    if sim_type == "stella":
        t, z, phi_vs_t, beta = extract_data_from_ncdf_with_xarray(outnc_longname, "t", 'zed', 'phi_vs_t', 'beta')
        phi_vs_t = np.array(phi_vs_t[:,0,:,0,0,:])
    elif sim_type == "gs2":
        t, z, phi_vs_t, beta = extract_data_from_ncdf_with_xarray(outnc_longname, "t", 'theta', 'phi_t', 'beta')
        phi_vs_t = np.array(phi_vs_t[:,0,0,:,:])
    zval= float(z[int(len(z)*0.5)])
    phi_t = phi_vs_t[:,int(len(z)*0.5),0]
    ([amp0_opt, phase_opt, freq_opt, gamma_opt, offset_opt],
           [amp0_err, phase_err, freq_err, gamma_err, offset_err]) = fit_ksaw(t, phi_t, **kwargs)

    return ([amp0_opt, phase_opt, freq_opt, gamma_opt, offset_opt],
           [amp0_err, phase_err, freq_err, gamma_err, offset_err])

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
    field_vs_t_list = []

    for idx, outnc_longname in enumerate(outnc_longnames):
        if sim_types[idx] == "stella":
            if field_name == "phi":
                field_key = "phi_vs_t"
            elif field_name == "apar":
                field_key = "apar_vs_t"
            elif field_name == "bpar":
                field_key = "bpar_vs_t"
            print("stella")
            view_ncdf_variables_with_xarray(outnc_longname)
            t, z, field_vs_t, beta = extract_data_from_ncdf_with_xarray(outnc_longname, "t", 'zed', field_key, 'beta')
            if field_name == "apar":
                field_vs_t = np.array(field_vs_t[:,0,:,0,0,:])
            else:
                field_vs_t = np.array(field_vs_t[:,0,:,0,0,:])
            # shape is now (time, zed, ri)
        elif sim_types[idx] == "gs2":
            if field_name == "phi":
                field_key = "phi_t"
            elif field_name == "apar":
                field_key = "apar_t"
            elif field_name == "bpar":
                field_key = "bpar_t"
            print("gs2")
            view_ncdf_variables_with_xarray(outnc_longname)
            t, z, field_vs_t, beta = extract_data_from_ncdf_with_xarray(outnc_longname, "t", 'theta', field_key, 'beta')
            field_vs_t = np.array(field_vs_t[:,0,0,:,:])
            if field_name == "apar":
                field_vs_t = field_vs_t/2
            # shape is now (time, zed, ri)
        else:
            print("sim type not recognised! Aborting")
            sys.exit()
        t_list.append(t); z_list.append(z); field_vs_t_list.append(field_vs_t);
    for idx, outnc_longname in enumerate(outnc_longnames):
        t = t_list[idx]; z = z_list[idx]; field_vs_t = field_vs_t_list[idx];
        zval= float(z[int(len(z)*0.5)])
        if abs(zval) > 1E-5:
            print("Error! zval!=0 . zval = ", zval)
            sys.exit()
        if field_name == "apar":
            field_t = field_vs_t[:,int(len(z)*0.25),0]
        else:
            field_t = field_vs_t[:,int(len(z)*0.5),0]
        if idx == 0:
            field_t_orig = field_t
            z_orig = zval
            t_orig = t
            field_t_normalising_factor = np.max(abs(field_t_orig))
        ax1.plot(t, field_t, label=(labels[idx]), ls=linestyles[idx], lw=my_linewidth)


        if len(t_orig) == len(t):
            if np.max(abs(t - t_orig)) < 1E-5:
                if idx != 0:    # Don't plot the first, as this is what we're comparing to
                    ax2.plot(t, (abs((field_t - field_t_orig)/field_t_normalising_factor)*100),
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
    gs2_implicit_dt4Em2 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar1/input_nz36_dt4E-2.out.nc"
    gs2_implicit_dt4Em3 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar1/input_nz36_dt4E-3.out.nc"
    gs2_implicit_dt4Em2_bd0 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar1/input_nz36_dt4E-2_bakdif0.out.nc"
    gs2_implicit_dt4Em2_bd0_fe05 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar1/input_nz36_dt4E-2_bakdif0_fexpr0.5.out.nc"
    gs2_implicit_dt4Em2_nz72 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar1/input_nz72_dt4E-2.out.nc"
    # stella_implicit_dt5Em2 = "sims/for_thesis_stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_implicit_nz36_dt5E-2.out.nc"
    stella_explicit_src_h_dt4Em4 = "sims/stella_comparing_implex_gbar_h_src/input_explicit_src_h.out.nc"
    stella_explicit_src_h_dt4Em4_zupwexpl01 = "sims/stella_comparing_implex_gbar_h_src/input_explicit_src_h_zupwexpl0.1.out.nc"
    stella_implicit_src_h_dt4Em4 = "sims/stella_comparing_implex_gbar_h_src/input_implicit_src_h.out.nc"
    stella_implicit_src_h_dt4Em2 = "sims/stella_comparing_implex_gbar_h_src/input_implicit_src_h_dt4E-2.out.nc"
    stella_implicit_src_h_exp = "sims/stella_comparing_implex_gbar_h_src/input_implicit_src_h_experimental.out.nc"
    stella_implicit_src_h_dt4Em2_tupw002 = "sims/stella_comparing_implex_gbar_h_src/input_implicit_src_h_dt4E-2_tupw002.out.nc"
    stella_implicit_src_h_dt4Em2_zupw002 = "sims/stella_comparing_implex_gbar_h_src/input_implicit_src_h_dt4E-2_zupw002.out.nc"

    # sim_types = [
    #             "gs2", "gs2", "gs2", "gs2", "gs2",
    #             "stella", "stella", "stella", #"stella",
    #             ]
    # for idx, outnc_longname in enumerate([gs2_implicit_dt4Em2,
    #             gs2_implicit_dt4Em3,
    #             gs2_implicit_dt4Em2_bd0,
    #             gs2_implicit_dt4Em2_bd0_fe05,
    #             gs2_implicit_dt4Em2_nz72,
    #             stella_explicit_src_h_dt4Em4,
    #             stella_implicit_src_h_dt4Em4,
    #             stella_implicit_src_h_dt4Em2,
    #             #stella_implicit_src_h_dt4Em2_zupw002
    #             ]):
    #     #######################
    #     ([amp0_opt, phase_opt, freq_opt, gamma_opt, offset_opt],
    #            [amp0_err, phase_err, freq_err, gamma_err, offset_err]) = fit_frequency_and_damping_rate(outnc_longname, sim_types[idx], make_plot=False)
    #     print("outnc_longname, gamma_opt = ", outnc_longname, gamma_opt)

    compare_sims_for_thesis("phi",
            [# gs2_implicit_dt25Em5,
             gs2_implicit_dt4Em2,
             # gs2_implicit_dt4Em3,
            # gs2_implicit_dt4Em2_bd0,
             # gs2_implicit_dt4Em2_bd0_fe05,
             # gs2_implicit_dt4Em2_nz72,
             stella_explicit_src_h_dt4Em4,
             stella_explicit_src_h_dt4Em4_zupwexpl01,
             stella_implicit_src_h_dt4Em2,
             #stella_implicit_src_h_exp,
             stella_implicit_src_h_dt4Em2_tupw002
             #stella_implicit_dt25Em5
             ],
            ["gs2",
             # "gs2",
             # "gs2",
             "stella", "stella", "stella", "stella"
             ],
            [
            # "gs2(dt=2.5E-4)",
            "gs2(dt=4E-2)",
            # "gs2(dt=4E-3)",
            #"gs2(dt=4E-2) (bakdif=0)",
            # "gs2(dt=4E-2) (bakdif=0, fexpr=0.5)",
            #"gs2(dt=4E-2), nz=72",
            "stella (explicit, src. h)",
            "stella (explicit, src. h, zupwexpl=0.1)",
            "stella (impl, src. h, dt=4E-2, zupw=tupw=0)",#, "stella (explicit, centered dg/dz)",
            #"stella (exp.)",
            "stella (impl, src. h, dt=4E-2, tupw=0.02)"#, "stella (explicit, centered dg/dz)",
            # "stella (dt=2.5E-4)",
            ],
             save_name="phi_t_fapar1_fbpar1_src_h.eps")

    compare_sims_for_thesis("bpar",
            [
             gs2_implicit_dt4Em2,
             stella_implicit_src_h_dt4Em2,
             ],
            ["gs2",
             "stella",
             ],
            [
            "gs2(dt=4E-2)",
            "stella (impl, src. h, dt=4E-2, zupw=tupw=0)",#, "stella (explicit, centered dg/dz)",
            ],
             save_name="bpar_t_fapar1_fbpar1_src_h.eps")

    compare_sims_for_thesis("apar",
            [# gs2_implicit_dt25Em5,
             gs2_implicit_dt4Em2,
             stella_implicit_src_h_dt4Em2,
             ],
            ["gs2",
             "stella"
             ],
            [
            "gs2(dt=4E-2)",
            "stella (impl, src. h, dt=4E-2, zupw=tupw=0)",#, "stella (explicit, centered dg/dz)",
            ],
             save_name="apar_t_fapar1_fbpar1_src_h.eps")

def benchamark_stella_src_h_stella_vs_gs2_fapar1_fbpar0():
    """ """
    gs2_implicit_dt5Em2 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar0/input_nz36_dt5E-2.out.nc"
    #stella_implicit_dt5Em2 = "sims/for_thesis_stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_implicit_nz36_dt5E-2.out.nc"
    stella_explicit_src_h_dt4Em4 = "sims/stella_comparing_implex_gbar_h_src/input_explicit_src_h_fapar1_fbpar0.out.nc"
    stella_explicit_src_h_dt4Em4_ = "sims/stella_comparing_implex_gbar_h_src/input_explicit_src_h_fapar1_fbpar0.out.nc"
    stella_explicit_src_h_dt4Em4_zupwexpl01 = "sims/stella_comparing_implex_gbar_h_src/input_explicit_src_h_fapar1_fbpar0_zupwexpl_0.1.out.nc"
    # stella_implicit_src_h_dt4Em4 = "sims/stella_comparing_implex_gbar_h_src/input_implicit_src_h.out.nc"
    stella_implicit_src_h_dt4Em2 = "sims/stella_comparing_implex_gbar_h_src/input_implicit_src_h_fapar1_fbpar0_dt4E-2.out.nc"


    compare_sims_for_thesis("phi",
            [# gs2_implicit_dt25Em5,
             gs2_implicit_dt5Em2,
             # stella_implicit_dt5Em2,
             stella_explicit_src_h_dt4Em4,
             #stella_implicit_src_h_dt4Em4,
             stella_implicit_src_h_dt4Em2,
             stella_explicit_src_h_dt4Em4_zupwexpl01
             #stella_implicit_dt25Em5
             ],
            ["gs2",
             "stella", "stella", "stella", # "stella"
             ],
            [
            # "gs2(dt=2.5E-4)",
            "gs2(dt=5E-4)",
            # "stella (implicit, dt=4E-2)",
            "stella (explicit, src. h)",
            "stella (explicit, src. h, zupwexpl=0.1)",
            #"stella (impl, src. h, dt=4E-4)",
            "stella (impl, src. h, dt=4E-2)"#, "stella (explicit, centered dg/dz)",

            # "stella (dt=2.5E-4)",
            ],
             save_name="for_thesis_phi_t_fapar1_fbpar0_src_h.eps")

    return

def benchmark_stella_src_h_stella_vs_gs2_fapar1_fbpar1_for_thesis():
    """
    Make a plot showing phi(z), apar(z), bpar(z), phi(t) for the fiducial stella implicit,
    stella explicit and GS2 sims.
    Then fit phi(t) to the KSAW and report Omega for the following sims:
     GS2
     - Fiducial GS2: gs2_implicit_dt4Em2
     - bakdif=0, fexpr=0.5
     - increased nzed
     - increased v-res
     - reduced dt (4E-4)
     stella, implicit
     - Fiducial: stella_implicit_src_h_dt4Em2
     - ut=0.02
     - uz=0.02 (unstable)
     - increased v-res
     - increased nzed
     - reduced dt (4E-4)
     stella, explicit
     - Fiducial: stella_explicit_src_h_dt4Em4_zupwexpl01
     - zupwexpl = 1
     - zupwexpl = 0
     - zupwexpl = 0.5
      """

    def make_fields_zt_plots_for_thesis():
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
        my_xticklength = 8
        my_xtickwidth = 2
        my_yticklength = 8
        my_ytickwidth = 2
        my_xticklength_minor = 6
        my_xtickwidth_minor = 1
        my_yticklength_minor = 6
        my_ytickwidth_minor = 1

        default_cols = plt.get_cmap("tab10")

        left = 0.17
        right = 0.985
        top = 0.985
        bottom = 0.08
        vspace=0.1
        hspace=0.15
        height = (top - bottom - vspace)/2
        width = (right - left - hspace)/2
        col2_left = left + width + hspace
        my_linewidth = 3

        fig = plt.figure(figsize=(12,12))
        ax1 = fig.add_axes((left, bottom + height + vspace, width, height))
        ax2 = fig.add_axes((col2_left, bottom + height + vspace, width, height))
        ax3 = fig.add_axes((left, bottom, width, height))
        ax4 = fig.add_axes((col2_left, bottom, width, height))

        outnc_longnames = fiducial_sims
        sim_types = fiducial_sim_types
        labels = fiducial_sim_labels
        save_name="for_thesis_phi_t_fapar1_fbpar1_src_h.eps"

        t_list = []
        z_list = []
        phi_vs_t_list = []
        apar_vs_t_list = []
        bpar_vs_t_list = []

        for idx, outnc_longname in enumerate(outnc_longnames):
            print("outnc_longname = ", outnc_longname)
            if sim_types[idx] == "stella":
                # view_ncdf_variables_with_xarray(outnc_longname)
                t, z, phi_vs_t, apar_vs_t, bpar_vs_t, beta = extract_data_from_ncdf_with_xarray(outnc_longname,
                    "t", 'zed', "phi_vs_t", "apar_vs_t", "bpar_vs_t", 'beta')
                phi_vs_t = phi_vs_t[:,0,:,0,0,:]
                apar_vs_t = apar_vs_t[:,0,:,0,0,:]
                bpar_vs_t = bpar_vs_t[:,0,:,0,0,:]
                # shape is now (time, zed, ri)
            elif sim_types[idx] == "gs2":
                # view_ncdf_variables_with_xarray(outnc_longname)
                t, z, phi_vs_t, apar_vs_t, bpar_vs_t, beta = extract_data_from_ncdf_with_xarray(outnc_longname, "t", 'theta', "phi_t",
                    "apar_t", "bpar_t", 'beta')
                phi_vs_t = phi_vs_t[:,0,0,:,:]
                apar_vs_t = apar_vs_t[:,0,0,:,:]/2  # normalise to stella normalisation
                bpar_vs_t = bpar_vs_t[:,0,0,:,:]
                # shape is now (time, zed, ri)
            else:
                print("sim type not recognised! Aborting")
                sys.exit()
            t_list.append(t); z_list.append(z);
            phi_vs_t_list.append(phi_vs_t) ;  apar_vs_t_list.append(apar_vs_t)
            bpar_vs_t_list.append(bpar_vs_t)

        ## Plot phi, apar, bpar at final timestep
        for idx, outnc_longname in enumerate(outnc_longnames):
            t = t_list[idx]; z = z_list[idx]; phi_vs_t = phi_vs_t_list[idx];
            apar_vs_t = apar_vs_t_list[idx]; bpar_vs_t = bpar_vs_t_list[idx]
            print("phi_vs_t.shape = ", apar_vs_t.shape)
            print("apar_vs_t.shape = ", apar_vs_t.shape)
            zval= float(z[int(len(z)*0.5)])
            if abs(zval) > 1E-5:
                print("Error! zval!=0 . zval = ", zval)
                sys.exit()
            apar_t = apar_vs_t[:,int(len(z)*0.25),0]
            phi_t = phi_vs_t[:,int(len(z)*0.5),0]
            # bpar_t = bpar_vs_t[:,int(len(z)*0.5),0]
            # print("phi_t = ", phi_t)
            # print("phi_vs_t[-1,:,0] = ", phi_vs_t[-1,:,0])
            ax1.plot(z/np.pi, phi_vs_t[-1,:,0], label=(labels[idx]), ls=linestyles[idx], lw=my_linewidth)
            ax2.plot(z/np.pi, apar_vs_t[-1,:,0], label=(labels[idx]), ls=linestyles[idx], lw=my_linewidth)
            ax3.plot(z/np.pi, bpar_vs_t[-1,:,0], label=(labels[idx]), ls=linestyles[idx], lw=my_linewidth)
            ax4.plot(t, phi_t, label=(labels[idx]), ls=linestyles[idx], lw=my_linewidth)
            #ax4.plot(t, apar_t, label=(labels[idx]), ls=linestyles[idx], lw=my_linewidth)

        for ax in [ax1, ax2, ax3]:
            ax.set_xlim((-1, 1))
            ax.plot([-1, 1], [0,0], c="black", lw=1., zorder=-10)
        ax4.set_xlim((0,40))
        ax4.plot([0,40], [0,0], c="black", lw=1., zorder=-10)
        # if field_name == "phi":
        #     field_ylabel_ax1 = r"$\tilde{\phi}_k$"
        #     field_ylabel_ax2 = r"$\tilde{\phi}_k$ error (%)"
        # if field_name == "apar":
        #     field_ylabel_ax1 = r"$\tilde{A}_{\parallel, k}$"
        #     field_ylabel_ax2 = r"$\tilde{A}_{\parallel, k}$ error (%)"
        # if field_name == "bpar":
        #     field_ylabel_ax1 = r"$\tilde{B}_{\parallel, k}$"
        #     field_ylabel_ax2 = r"$\tilde{B}_{\parallel, k}$ error (%)"

        # ax1.set_xlim((-1, 1))
        # ax2.set_xlim((-1, 1))
        # ax1.set_xlim((0,20))
        # ax1.set_ylim((-6, 6))
        # if ax1_yticks is not None:
        #     ax1.set_yticks(ax1_yticks)
        #     ax1.set_yticklabels(ax1_yticklabels, fontsize=y_ticklabelfontsize)
        # if ax2_yticks is not None:
        #     ax2.set_yticks(ax2_yticks)
        #     ax2.set_yticklabels(ax2_yticklabels, fontsize=y_ticklabelfontsize)
        for ax in [ax1, ax2, ax3, ax4]:
            ax.tick_params("x", length=my_xticklength, width=my_xtickwidth, direction="in")
            ax.tick_params("y", length=my_yticklength, width=my_ytickwidth, direction="in")
            ax.tick_params("x", which="minor", length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in")
            ax.tick_params("y", which="minor", length=my_yticklength_minor, width=my_ytickwidth_minor, direction="in")
        for ax in [ax1, ax2, ax3]:
            ax.set_xticks([-1, 0, 1])
            ax.set_xticklabels([r"$-\pi$", r"$0$", r"$\pi$"], fontsize=x_ticklabelfontsize)
            ax.set_xticks([-0.5, 0.5], minor=True)
            ax.set_xticklabels([], minor=True)

        ax1.set_yticks([-0.2, 0, 0.2])
        ax1.set_yticklabels([r"$-0.2$", r"$0$", r"$0.2$"], fontsize=y_ticklabelfontsize)
        ax1.set_yticks([-0.3, -0.1, 0.1, 0.3], minor=True)
        ax1.set_yticklabels([], minor=True)

        ax2.set_yticks([-1, 0, 1])
        ax2.set_yticklabels([r"$-0.2$", r"$0$", r"$0.2$"], fontsize=y_ticklabelfontsize)
        ax2.set_yticks([-0.5, 0.5], minor=True)
        ax2.set_yticklabels([], minor=True)

        ax3.set_yticks([-0.04, 0, 0.04])
        ax3.set_yticklabels([r"$-0.04$", r"$0$", r"$0.04$"], fontsize=y_ticklabelfontsize)
        ax3.set_yticks([-0.06, -0.02, 0.02, 0.06], minor=True)
        ax3.set_yticklabels([], minor=True)

        ax4.set_yticks([-4, -2, 0, 2])
        ax4.set_yticklabels([r"$-4$", r"$-2$", r"$0$", r"$2$"], fontsize=y_ticklabelfontsize)
        ax4.set_yticks([-3, -1, 1, 3], minor=True)
        ax4.set_yticklabels([], minor=True)

        ax2.legend(loc="best", fontsize=legend_fontsize)
        ax1.set_ylabel(r"$\tilde{\phi}_k(\tilde{t}=40)$", fontsize=ylabel_fontsize)
        ax2.set_ylabel(r"$\tilde{A}_{\parallel,k}(\tilde{t}=40)$", fontsize=ylabel_fontsize, labelpad=-10)
        ax3.set_ylabel(r"$\tilde{B}_{\parallel,k}(\tilde{t}=40)$", fontsize=ylabel_fontsize, labelpad=0)
        ax4.set_ylabel(r"$\tilde{\phi}_k(z=0)$", fontsize=ylabel_fontsize)
        ax4.set_xlabel(r"$\tilde{t}$", fontsize=xlabel_fontsize)
        for ax in [ax1, ax2, ax3]:
            ax.set_xlabel(r"$z$", fontsize=xlabel_fontsize)

        # plt.show()
        plt.savefig(save_name)
        plt.close()


        return



    ## GS2
    gs2_implicit_dt4Em2 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar1/input_nz36_dt4E-2.out.nc"
    gs2_implicit_dt4Em2_bd0_fe05 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar1/input_nz36_dt4E-2_bakdif0_fexpr0.5.out.nc"
    gs2_implicit_dt4Em2_nz72 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar1/input_nz72_dt4E-2.out.nc"
    gs2_implicit_dt4Em4 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar1/input_nz36_dt4E-4.out.nc"
    gs2_implicit_dt4Em2_higher_vres = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar1/input_nz36_dt4E-2_higher_vres.out.nc"
    gs2_implicit_dt4Em2_quad_vres = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar1/input_nz36_dt4E-2_quad_vres.out.nc"
    # stella implicit
    stella_implicit_src_h_dt4Em2 = "sims/stella_src_h_for_thesis/input_implicit_src_h_dt4E-2.out.nc"
    stella_implicit_src_h_dt4Em4 = "sims/stella_src_h_for_thesis/input_implicit_src_h.out.nc"
    stella_implicit_src_h_dt4Em2_tupw002 = "sims/stella_src_h_for_thesis/input_implicit_src_h_dt4E-2_tupw002.out.nc"
    stella_implicit_src_h_dt4Em2_zupw002 = "sims/stella_src_h_for_thesis/input_implicit_src_h_dt4E-2_zupw002.out.nc"
    stella_implicit_src_h_dt4Em2_nz72 = "sims/stella_src_h_for_thesis/input_implicit_src_h_dt4E-2_nzed72.out.nc"
    stella_implicit_src_h_dt4Em2_higher_vres = "sims/stella_src_h_for_thesis/input_implicit_src_h_dt4E-2_higher_vres.out.nc"
    # stella explicit
    stella_explicit_src_h_dt4Em4_zupwexpl01 = "sims/stella_src_h_for_thesis/input_explicit_src_h_zupwexpl0.1.out.nc"
    stella_explicit_src_h_dt4Em4_zupwexpl1 = "sims/stella_src_h_for_thesis/input_explicit_src_h.out.nc"
    stella_explicit_src_h_dt4Em4_zupwexpl0 = "sims/stella_src_h_for_thesis/input_explicit_src_h_zupwexpl0.out.nc"
    stella_explicit_src_h_dt4Em4_zupwexpl05 = "sims/stella_src_h_for_thesis/input_explicit_src_h_zupwexpl0.5.out.nc"

    fiducial_sims = [gs2_implicit_dt4Em2, stella_implicit_src_h_dt4Em2,
                    # stella_explicit_src_h_dt4Em4_zupwexpl01
                    ]
    fiducial_sim_types = ["gs2", "stella",
                        # "stella"
                        ]
    fiducial_sim_labels = ["GS2 (fid.)", "stella impl. (fid.)",
                        # "stella explicit (fid.)"
                            ]
    make_fields_zt_plots_for_thesis()

    # for idx, outnc_longname in enumerate([gs2_implicit_dt4Em2,
    #         gs2_implicit_dt4Em2_bd0_fe05,
    #         gs2_implicit_dt4Em2_nz72,
    #         gs2_implicit_dt4Em4,
    #         gs2_implicit_dt4Em2_higher_vres,
    #         gs2_implicit_dt4Em2_quad_vres,
    #              ]):
    #     #######################
    #     ([amp0_opt, phase_opt, freq_opt, gamma_opt, offset_opt],
    #      [amp0_err, phase_err, freq_err, gamma_err, offset_err]) = fit_frequency_and_damping_rate(outnc_longname,
    #                         "gs2", make_plot=True)
    #     print("outnc_longname  = ", outnc_longname)
    #     print("(freq, err, gamma, err) = ",
    #         ("(" + str(freq_opt) + ", " + str(freq_err) + ", " + str(gamma_opt) + ", " + str(gamma_err) + ")"))

    # for idx, outnc_longname in enumerate([stella_implicit_src_h_dt4Em2,
    #         stella_implicit_src_h_dt4Em4,
    #         stella_implicit_src_h_dt4Em2_tupw002,
    #         stella_implicit_src_h_dt4Em2_nz72,
    #         stella_implicit_src_h_dt4Em2_higher_vres
    #              ]):
    #     #######################
    #     ([amp0_opt, phase_opt, freq_opt, gamma_opt, offset_opt],
    #      [amp0_err, phase_err, freq_err, gamma_err, offset_err]) = fit_frequency_and_damping_rate(outnc_longname,
    #                         "stella", make_plot=True)
    #     print("outnc_longname  = ", outnc_longname)
    #     print("(freq, err, gamma, err) = ",
    #         ("(" + str(freq_opt) + ", " + str(freq_err) + ", " + str(gamma_opt) + ", " + str(gamma_err) + ")"))

    for idx, outnc_longname in enumerate([stella_explicit_src_h_dt4Em4_zupwexpl01,
        stella_explicit_src_h_dt4Em4_zupwexpl0,
        stella_explicit_src_h_dt4Em4_zupwexpl05,
        stella_explicit_src_h_dt4Em4_zupwexpl1,
                 ]):
        #######################
        ([amp0_opt, phase_opt, freq_opt, gamma_opt, offset_opt],
         [amp0_err, phase_err, freq_err, gamma_err, offset_err]) = fit_frequency_and_damping_rate(outnc_longname,
                            "stella", make_plot=True)
        print("outnc_longname  = ", outnc_longname)
        print("(freq, err, gamma, err) = ",
            ("(" + str(freq_opt) + ", " + str(freq_err) + ", " + str(gamma_opt) + ", " + str(gamma_err) + ")"))

    return

def benchamark_stella_src_h_stella_vs_gs2_fapar1_fbpar0_for_thesis():
    """ """
    gs2_implicit_dt5Em2 = "sims/for_thesis_gs2_ky1_beta1_zero_gradients_fapar1_fbpar0/input_nz36_dt5E-2.out.nc"
    #stella_implicit_dt5Em2 = "sims/for_thesis_stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_implicit_nz36_dt5E-2.out.nc"
    stella_explicit_src_h_dt4Em4 = "sims/stella_comparing_implex_gbar_h_src/input_explicit_src_h_fapar1_fbpar0.out.nc"
    stella_explicit_src_h_dt4Em4_ = "sims/stella_comparing_implex_gbar_h_src/input_explicit_src_h_fapar1_fbpar0.out.nc"
    stella_explicit_src_h_dt4Em4_zupwexpl01 = "sims/stella_comparing_implex_gbar_h_src/input_explicit_src_h_fapar1_fbpar0_zupwexpl_0.1.out.nc"
    # stella_implicit_src_h_dt4Em4 = "sims/stella_comparing_implex_gbar_h_src/input_implicit_src_h.out.nc"
    stella_implicit_src_h_dt4Em2 = "sims/stella_comparing_implex_gbar_h_src/input_implicit_src_h_fapar1_fbpar0_dt4E-2.out.nc"


    compare_sims_for_thesis("phi",
            [# gs2_implicit_dt25Em5,
             gs2_implicit_dt5Em2,
             # stella_implicit_dt5Em2,
             stella_explicit_src_h_dt4Em4,
             #stella_implicit_src_h_dt4Em4,
             stella_implicit_src_h_dt4Em2,
             stella_explicit_src_h_dt4Em4_zupwexpl01
             #stella_implicit_dt25Em5
             ],
            ["gs2",
             "stella", "stella", "stella", # "stella"
             ],
            [
            # "gs2(dt=2.5E-4)",
            "gs2(dt=5E-4)",
            # "stella (implicit, dt=4E-2)",
            "stella (explicit, src. h)",
            "stella (explicit, src. h, zupwexpl=0.1)",
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
    # benchmark_stella_src_h_stella_vs_gs2_fapar1_fbpar1()
    # benchamark_stella_src_h_stella_vs_gs2_fapar1_fbpar0()
    benchmark_stella_src_h_stella_vs_gs2_fapar1_fbpar1_for_thesis()
    # benchamark_stella_src_h_stella_vs_gs2_fapar1_fbpar0_for_thesis()
