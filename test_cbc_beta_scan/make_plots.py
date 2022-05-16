""" """

import sys
sys.path.append("../postprocessing_tools")
from plotting_helper import make_beta_scan_plots, make_comparison_plots, plot_gmvus, plot_gzvs
from plotting_helper import make_comparison_plots_for_poster
from helper_ncdf import view_ncdf_variables, extract_data_from_ncdf
from extract_sim_data import find_bpar_phi_ratio, find_apar_phi_ratio
import matplotlib.pyplot as plt
import numpy as np

IMAGE_DIR = "./images/"
### Give the names of all the sims here - avoids needing to type them out
### in the methods.

## fapar = 0
# stella
sim_st_b00001_fapar0 = "stella_beta0.00001_fapar0/input"
sim_st_b001_fapar0 = "stella_beta0.001_fapar0/input"
sim_st_b002_fapar0 = "stella_beta0.002_fapar0/input"
sim_st_b01_fapar0  = "stella_beta0.010_fapar0/input"
# gs2
sim_gs2_b00001_fapar0 = "gs2_beta_scan_fapar0/_0.00001"
sim_gs2_b001_fapar0 = "gs2_beta_scan_fapar0/_0.0010"
sim_gs2_b002_fapar0 = "gs2_beta_scan_fapar0/_0.0020"
sim_gs2_b01_fapar0 = "gs2_beta_scan_fapar0/_0.0100"

## fbpar = 0
# stella
sim_st_b00001_fbpar0 = "stella_beta0.00001_fbpar0/input"
sim_st_b00005_fbpar0 = "stella_beta0.00005_fbpar0/input"
sim_st_b0001_fbpar0 = "stella_beta0.0001_fbpar0/input"
sim_st_b0003_fbpar0 = "stella_beta0.0003_fbpar0/input"
sim_st_b0006_fbpar0 = "stella_beta0.0006_fbpar0/input"
sim_st_b0015_fbpar0 = "stella_beta0.0015_fbpar0/input"
sim_st_b002_fbpar0 = "stella_beta0.002_fbpar0/input"
sim_st_b003_fbpar0 = "stella_beta0.003_fbpar0/input"
sim_st_b004_fbpar0 = "stella_beta0.004_fbpar0/input"
sim_st_b005_fbpar0 = "stella_beta0.005_fbpar0/input"
sim_st_b01_fbpar0 = "stella_beta0.010_fbpar0/input"
#######
sim_st_b001_fbpar0 = "stella_beta0.001_fbpar0/input"
sim_st_b001_fbpar0_np2 = "stella_beta0.001_fbpar0/input_np2"
sim_st_b001_fbpar0_np6 = "stella_beta0.001_fbpar0/input_np6"
sim_st_b001_fbpar0_nzed32 = "stella_beta0.001_fbpar0/input_nzed32"
sim_st_b001_fbpar0_nzed128 = "stella_beta0.001_fbpar0/input_nzed128"
sim_st_b001_fbpar0_no_drive = "stella_beta0.001_fbpar0/input_zero_drive"
sim_st_b001_fbpar0_no_drive_no_stream_no_mirror = "stella_beta0.001_fbpar0/input_zero_drive_zero_streaming_zero_mirror"
sim_st_b001_fbpar0_no_stream_no_mirror = "stella_beta0.001_fbpar0/input_zero_streaming_zero_mirror"
sim_st_b001_fbpar0_no_drive_0_upwind = "stella_beta0.001_fbpar0/input_zero_drive_no_upwinding"
sim_st_b001_fbpar0_no_drive_01_upwind = "stella_beta0.001_fbpar0/input_zero_drive_more_upwinding"
sim_st_b001_fbpar0_no_drive_02_upwind = "stella_beta0.001_fbpar0/input_zero_drive_0.2_upwinding"
sim_st_b001_fbpar0_no_drive_mid_vres = "stella_beta0.001_fbpar0/input_zero_drive_mid_vres"
sim_st_b001_fbpar0_no_drive_higher_vres = "stella_beta0.001_fbpar0/input_zero_drive_higher_vres"
# NB _no_drive_no_drifts has diamagnetic and curvature drifts set to zero
sim_st_b001_fbpar0_no_drive_no_drifts = "stella_beta0.001_fbpar0/input_zero_drive_zero_drifts"
sim_st_b001_fbpar0_no_drive_lower_dt = "stella_beta0.001_fbpar0/input_zero_drive_lower_dt"
# NB has diamagnetic but no curvature drifts
sim_st_b001_fbpar0_no_mag_drift = "stella_beta0.001_fbpar0/input_no_drifts"
sim_st_b001_fbpar0_no_mirror = "stella_beta0.001_fbpar0/input_no_mirror"
sim_st_b001_fbpar0_no_streaming = "stella_beta0.001_fbpar0/input_no_streaming"
sim_st_b001_fbpar0_no_zed_upwind = "stella_beta0.001_fbpar0/input_no_zed_upwinding" ## Doesn't actually do anything
sim_st_b001_fbpar0_centered_dgdz = "stella_beta0.001_fbpar0/input_centered_dgdz"
sim_st_b001_fbpar0_centered_dgdvpa = "stella_beta0.001_fbpar0/input_centered_dgdz_and_dgdvpa"
sim_st_b001_fbpar0_centered_dgdvpa_dgdvpa = "stella_beta0.001_fbpar0/input_centered_dgdz_and_dgdvpa"
sim_st_b001_fbpar0_centered_dgdvpa_dgdvpa_numapar = "stella_beta0.001_fbpar0/input_centered_dgdz_and_dgdvpa_num_apar_fac"
sim_st_b001_fbpar0_equal_masses = "stella_beta0.001_fbpar0_equal_masses/input"
sim_st_b001_fbpar0_flipflop = "stella_beta0.001_fbpar0/input_np5_flipflop"
###################
# gs2
sim_gs2_b00001_fbpar0 = "gs2_beta_scan_fbpar0/_0.00001"
sim_gs2_b001_fbpar0 = "gs2_beta_scan_fbpar0/_0.0010"
sim_gs2_b001_fbpar0_equal_masses = "gs2_beta_scan_fbpar0/_0.0010_me1"
sim_gs2_b002_fbpar0 = "gs2_beta_scan_fbpar0/_0.0020"
sim_gs2_b003_fbpar0 = "gs2_beta_scan_fbpar0/_0.0030"
sim_gs2_b004_fbpar0 = "gs2_beta_scan_fbpar0/_0.0040"
sim_gs2_b005_fbpar0 = "gs2_beta_scan_fbpar0/_0.0050"
sim_gs2_b01_fbpar0 = "gs2_beta_scan_fbpar0/_0.0100"
sim_gs2_b03 = "gs2_beta0.03/_0.0300"
sim_gs2_b03_hvr = "gs2_beta0.03/_0.0300_higher_vres"
sim_gs2_b01 = "gs2_beta0.01/_0.0100"
sim_gs2_b03_fbpar0 = "gs2_beta0.03_fbpar0/_0.0300"
sim_gs2_b03_fapar0 = "gs2_beta0.03_fapar0/_0.0300"

## fapar=1, fbpar=1
sim_st_b0 = "stella_beta0.000/input"
sim_st_b005 = "stella_beta0.005/input"
sim_st_b01 = "stella_beta0.010/input"
sim_st_b01_fbpar0_new = "stella_fapar1_fbpar0_beta_scan/beta_0.01000"
sim_st_b015 = "stella_beta0.015/input"
sim_st_b02 = "stella_beta0.020/input"
sim_st_b025 = "stella_beta0.025/input"
#sim_st_b03 = "stella_beta0.030/input"
sim_st_b03 = "stella_beta0.03/beta_0.03000"
sim_st_b03_fapar0 = "stella_beta0.03_fapar0/beta_0.03000"
sim_st_b03_fbpar0 = "stella_beta0.03_fbpar0/beta_0.03000"

### Local sims
sim_st_gfvmulo_nstep100 = "local_sims/stella_lr_gfvmulo_nstep100"
sim_st_gf_nstep100 = "local_sims/stella_lr_gf_nstep100"
sim_st_gfvmulo_nstep1000 = "local_sims/stella_lr_gfvmulo_beta1_nstep1000"
sim_st_gf_nstep1000 = "local_sims/stella_lr_gf_beta1_nstep1000"

## The pickled files summarising G2 beta scans
pickle_gs2 = "gs2_beta_scan/omega_values.pickle"
pickle_gs2_fbpar0 = "gs2_beta_scan_fbpar0/omega_values.pickle"

def analyse_fbpar0_beta0001_results():
    """Compare sims, all with fbpar=0, fapar=1, beta=0.001, for which
    we try turning on and off different knobs."""

    make_comparison_plots([
                           sim_st_b001_fbpar0,
                           sim_st_b001_fbpar0_no_drive,
                           sim_st_b001_fbpar0_no_mag_drift,
                           sim_st_b001_fbpar0_no_mirror,
                           sim_st_b001_fbpar0_equal_masses,
                           sim_st_b001_fbpar0_flipflop,
                           sim_st_b001_fbpar0_no_streaming,
                           sim_gs2_b001_fbpar0
                           ],
                          [
                           "stella",
                           "stella, zero gradients",
                           "stella, no magnetic drifts",
                           "stella, no mirror term",
                           "stella, m_e=1",
                           "stella, flip-flop",
                           "stella, no streaming",
                           "GS2"
                           ],
                          "images/termsoff_beta_0.001_fbpar0",
                          sim_types=[
                                     "stella",
                                     "stella",
                                     "stella",
                                     "stella",
                                     "stella",
                                     "stella",
                                     "stella",
                                    # "stella",
                                     "gs2"
                                     ],
                           plot_apar=True,
                           plot_format=".eps"
                           )
    make_comparison_plots([
                           sim_st_b001_fbpar0,
                            sim_st_b001_fbpar0_np2,
                            sim_st_b001_fbpar0_np6,
                            sim_gs2_b001_fbpar0
                           ],
                          [
                           "stella",
                           "stella, nperiod=2",
                           "stella, nperiod=6",
                           "GS2"
                           ],
                          "images/beta_0.001_fbpar0_nperiod_scan",
                          sim_types=[
                                     "stella",
                                     "stella",
                                     "stella",
                                     "gs2"
                                     ],
                           plot_apar=True,
                           plot_format=".eps"
                           )
    make_comparison_plots([
                           sim_st_b001_fbpar0,
                            sim_st_b001_fbpar0_nzed32,
                            sim_st_b001_fbpar0_nzed128,
                            sim_gs2_b001_fbpar0
                           ],
                          [
                           "stella",
                           "stella, nzed=32",
                           "stella, nzed=128",
                           "GS2"
                           ],
                          "images/beta_0.001_fbpar0_nzed_scan",
                          sim_types=[
                                     "stella",
                                     "stella",
                                     "stella",
                                     "gs2"
                                     ],
                           plot_apar=True,
                           plot_format=".png"
                           )
    make_comparison_plots([
                           sim_st_b001_fbpar0,
                            sim_st_b001_fbpar0_centered_dgdz,
                            sim_st_b001_fbpar0_centered_dgdvpa_dgdvpa,
                            #sim_st_b001_fbpar0_centered_dgdvpa_dgdvpa_numapar,
                            sim_gs2_b001_fbpar0
                           ],
                          [
                           "stella, dg/dz third-order-upwind",
                           "stella, dg/dz centered on proc0",
                           "stella, dg/dz centered on all procs",
                           #"stella, dg/dz centered, dg/dvpa centered, num.",
                           "GS2"
                           ],
                          "images/beta_0.001_fbpar0_dgdz_centering",
                          sim_types=[
                                     "stella",
                                     "stella",
                                     "stella",
                                     #"stella",
                                     "gs2"
                                     ],
                           plot_apar=True,
                           plot_format=".png"
                           )
    compare_omega_for_fbpar0_changing_streaming_and_drive()

    return

def analyse_results_for_poster():
    """Compare sims, all with fbpar=1, fapar=1, beta=0.01, for which
    we try turning on and off different knobs."""

    # make_comparison_plots_for_poster([
    #                        sim_st_b01,
    #                        sim_gs2_b01
    #                        ],
    #                       [
    #                        "stella",
    #                        "GS2"
    #                        ],
    #                       "images/beta_0.01_poster",
    #                       sim_types=[
    #                                  "stella",
    #                                  "gs2"
    #                                  ],
    #
    #                        )

    make_comparison_plots_for_poster([
                            sim_st_b03,
                            sim_gs2_b03,
                            #sim_gs2_b03_hvr,
                         ],
                        [
                         "stella",
                         "gs2",
                        # "gs2, hvr",
                         ],
                        IMAGE_DIR + "beta=0.03_poster",
                        sim_types=[
                                   "stella",
                                   "gs2",
                                   #"gs2",
                                   ],

                         )


    #compare_omega_for_fbpar0_changing_streaming_and_drive()

    return

def plot_fbpar0_beta0001_equal_masses():
    """Compare sims, all with fbpar=0, fapar=1, beta=0.001, for which
    we try turning on and off different knobs."""

    make_comparison_plots([
                           sim_st_b001_fbpar0,
                           sim_gs2_b001_fbpar0,
                           sim_st_b001_fbpar0_equal_masses,
                           sim_gs2_b001_fbpar0_equal_masses,
                           ],
                          [
                           "stella",
                           "GS2",
                           "stella, m_e=1",
                           "GS2, m_e=1",
                           ],
                          "./termsoff_beta_0.001_fbpar0_me1",
                          sim_types=[
                                     "stella",
                                     "gs2",
                                     "stella",
                                     "gs2",
                                     ],
                           plot_apar=True,
                           plot_format=".eps"
                           )
    return

def analyse_fbpar0_results():
    """Compare omega(t) and phi(z), apar(z)
    between GS2 and stella results """

    print("Hello world")
    beta_strs = [

                 ]
    stella_sim_longnames = [

                            ]
    gs2_sim_longnames = [

                            ]
    for beta_idx in range(0, len(beta_strs)):
        stella_sim_longname = stella_sim_longnames[beta_idx]
        gs2_sim_longname = gs2_sim_longnames[beta_idx]
        beta_str = beta_strs[beta_idx]
        make_comparison_plots([
                               stella_sim_longname,
                               gs2_sim_longname,
                               ],
                              [
                               "stella",
                               "GS2",
                               ],
                              "./beta_" + beta_str + "_fbpar0",
                              sim_types=[
                                         "stella",
                                         "gs2",
                                         ],
                               plot_apar=True,
                               )
    return

def analyse_fapar0_results():
    """Compare omega(t) and phi(z), apar(z)
    between GS2 and stella results """

    print("Hello world")
    beta_strs = [
                 "0.00001",
                 "0.001",
                 "0.002",
                 #"0.010",
                 ]
    stella_sim_longnames = [

                            ]
    gs2_sim_longnames = [

                            ]
    for beta_idx in range(0, len(beta_strs)):
        stella_sim_longname = stella_sim_longnames[beta_idx]
        gs2_sim_longname = gs2_sim_longnames[beta_idx]
        beta_str = beta_strs[beta_idx]
        make_comparison_plots([
                               stella_sim_longname,
                               gs2_sim_longname,
                               ],
                              [
                               "stella",
                               "GS2",
                               ],
                              "./beta_" + beta_str + "_fapar0",
                              sim_types=[
                                         "stella",
                                         "gs2",
                                         ],
                               plot_bpar=True,
                               plot_format=".png"
                               )
        # plot_gmvus("stella_beta0.001_fapar0/input.out.nc")
        print("stella_sim_longname = ", stella_sim_longname)
        print("gs2_sim_longname = ", gs2_sim_longname)
        gs2_bpar_ratio = find_bpar_phi_ratio(gs2_sim_longname, "gs2")
        stella_bpar_ratio = find_bpar_phi_ratio(stella_sim_longname, "stella")
        print("gs2_bpar_ratio = ", gs2_bpar_ratio)
        print("stella_bpar_ratio2 = ", stella_bpar_ratio)

    return

def analyse_fapar0_changing_vpares():
    """Compare omega(t) and phi(z), apar(z)
    between GS2 and stella results """

    print("Hello world")

    stella_sim1_longname = "stella_beta0.010_fapar0/input"
    stella_sim2_longname = "stella_beta0.010_fapar0_higher_vpa/input"

    gs2_sim1_longname = "gs2_beta_scan_fapar0/_0.0100"

    # make_comparison_plots([
    #                        stella_sim1_longname,
    #                        stella_sim2_longname,
    #                        gs2_sim1_longname
    #                        ],
    #                       [
    #                        "stella",
    #                        "stella, higher vpa",
    #                        "GS2",
    #                        ],
    #                       "./beta_0.01_fapar0_varying_vpares",
    #                       sim_types=[
    #                                  "stella",
    #                                  "stella",
    #                                  "gs2",
    #                                  ],
    #                        plot_bpar=True,
    #                        )
    #
    gs2_bpar_ratio = find_bpar_phi_ratio(gs2_sim1_longname, "gs2")
    stella_bpar_ratio1 = find_bpar_phi_ratio(stella_sim1_longname, "stella")
    stella_bpar_ratio2 = find_bpar_phi_ratio(stella_sim2_longname, "stella")
    print("gs2_bpar_ratio = ", gs2_bpar_ratio)
    print("stella_bpar_ratio1 = ", stella_bpar_ratio1)
    print("stella_bpar_ratio2 = ", stella_bpar_ratio2)
    # plot_gmvus("stella_beta0.001_fapar0/input.out.nc", which="gvpa")
    # plot_gmvus("stella_beta0.010_fapar0_higher_vpa/input.out.nc", which="gvpa")
    return

def plot_g():
    """ """
    stella_outnc_longname = "stella_beta0.001_fbpar0/input.out.nc"
    gs2_outnc_longname = "gs2_beta_scan_fbpar0/_0.0010.out.nc"
    #view_ncdf_variables(stella_outnc_longname)
    #view_ncdf_variables(gs2_outnc_longname)
    plot_gmvus(stella_outnc_longname)

    sys.exit()

    time, theta, gs2_energy, gs2_lambda = extract_data_from_ncdf(gs2_outnc_longname,
                                     't', 'theta', 'energy', 'lambda')
    print("len(time), len(theta), len(energy), len(lambda) = ", len(time), len(theta), len(gs2_energy), len(gs2_lambda))
    gs2_g = np.loadtxt("gs2_beta_scan_fbpar0/_0.0100.dist", skiprows=1)
    print("gs2_g.shape = ", gs2_g.shape)
    # GS2's g is in a set of (n(energy)*n(lambda)) x 8 blocks.
    # Each row is vpa, vpe, energy(ie,is), al(il), xpts(ie,is), ypts(il), real(gtmp(1)), real(gtmp(2))
    # The number of blocks is nstep/(nwrite*nwrite_mul)
    block_size = len(gs2_energy) * len(gs2_lambda)
    nblocks = len(gs2_g)/block_size
    final_step_g = gs2_g[-block_size:, :]
    print("len(final_step_g), block_size = ", len(final_step_g), block_size)

    ## Code to plot g for GS2
    # gvmus = gvmus[-1]   # spec, mu, vpa
    # fig = plt.figure()
    # ax1 = fig.add_subplot(211)
    # ax2 = fig.add_subplot(212)
    # counter=0
    #
    # for mu_idx in range(0, len(mu)):
    #     counter += 1
    #     g_ion_vpa = gvmus[0, mu_idx, :]
    #     g_electron_vpa = gvmus[1, mu_idx, :]
    #     ax1.plot(vpa, g_ion_vpa)
    #     ax2.plot(vpa, g_electron_vpa)
    #
    #     if counter == 5:
    #         plt.show()
    #         fig = plt.figure()
    #         ax1 = fig.add_subplot(211)
    #         ax2 = fig.add_subplot(212)
    #         counter=0
    #
    # plt.show()

    return

def plot_geometry():
    """ """
    stella_outnc_longname = "stella_beta0.001_fbpar0/input.out.nc"
    gs2_outnc_longname = "gs2_beta_scan_fbpar0/_0.0010.out.nc"

    z, gds2, gds21, gds22, bmag, gradpar = extract_data_from_ncdf(stella_outnc_longname,
                                    'zed', 'gds2', 'gds21', 'gds22', 'bmag', 'gradpar')

    # Code to compare geometry between stella and gs2
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax1.plot(z, gds2)
    ax1.plot(z, gds21)
    ax1.plot(z, gds22)
    ax2.plot(z, bmag)
    ax2.plot(z, gradpar)

    theta, gds2, gds21, gds22, bmag, gradpar = extract_data_from_ncdf(gs2_outnc_longname,
                                    'theta', 'gds2', 'gds21', 'gds22', 'bmag', 'gradpar')
    ax1.plot(theta, gds2, linestyle="-.")
    ax1.plot(theta, gds21, linestyle="-.")
    ax1.plot(theta, gds22, linestyle="-.")
    ax2.plot(theta, bmag, linestyle="-.")
    ax2.plot(theta, gradpar, linestyle="-.")


    plt.show()

    return

def plot_beta_scans():
    """ """
    stella_sim_longnames = [
                        sim_st_b0,
                        sim_st_b005,
                        sim_st_b01,
                        sim_st_b015,
                        sim_st_b02,
                        sim_st_b025,
                        sim_st_b03,
                            ]
    gs2_sim_longnames = [

                        ]
    stella_beta_vals = [
                        0.0,
                        0.005,
                        0.01,
                        0.015,
                        0.02,
                        0.025,
                        0.03
                        ]


    make_beta_scan_plots(stella_sim_longnames,
                         stella_beta_vals,
                          "images/test_cbc_beta_scan",
                         )

    # stella_sim_longnames = [
    #
    #                         ]
    # stella_beta_vals = [
    #                     0.0,
    #                     0.001,
    #                     0.002,
    #                     0.003,
    #                     0.004,
    #                     0.005,
    #                     ]
    #
    #
    #
    # make_beta_scan_plots(stella_sim_longnames,
    #                      stella_beta_vals,
    #                       "./test_cbc_beta_scan_fbpar0",
    #                      gs2_pickle=gs2_pickle
    #                      )
    return

def compare_vres_fbpar0():
    """ """
    stella_sim_longname = "stella_beta0.010_fbpar0/input"
    stella_sim_longname_higher_vpa = "stella_beta0.010_fbpar0_higher_vpa/input"
    stella_outnc_longname = stella_sim_longname + ".out.nc"
    stella_outnc_longname_higher_vpa = stella_sim_longname_higher_vpa + ".out.nc"
    make_comparison_plots([
                stella_sim_longname,
                stella_sim_longname_higher_vpa,
                ],
                [
                "nvgrid=36",
                "nvgrid=108",
                ],
                "./beta0.01_fbpar0_vpa_res_test",
                sim_types=[
                "stella",
                "stella",
                ],
                plot_apar=True,
                )
    return

def compare_omega_for_fbpar0_zero_drive():
    """ """

    make_comparison_plots([
                        sim_st_b001_fbpar0_no_drive,
                        sim_st_b001_fbpar0_no_drive_mid_vres,
                        #sim_st_b001_fbpar0_no_drive_higher_vres,
                        sim_st_b001_fbpar0_no_drive_lower_dt,
                        sim_gs2_b001_fbpar0,
                           ],
                          [
                           "stella, zero drive",
                           "stella, zero drive, mid vres",
                           #"stella, zero drive, higher vres",
                           "stella, zero drive, lower dt",
                           "GS2",
                           ],
                          IMAGE_DIR + "fbpar0_beta1e-3_zerodrive",
                          sim_types=[
                                     "stella",
                                     "stella",
                                     "stella",
                                     #"stella",
                                     "gs2",
                                     ],
                           plot_apar=True,
                           plot_bpar=False,
                           plot_format=".png"
                           )

def compare_omega_for_fbpar0_zero_drive_change_upwind():
    """ """

    make_comparison_plots([
                        sim_st_b001_fbpar0_no_drive,
                        sim_st_b001_fbpar0_no_drive_0_upwind,
                        sim_st_b001_fbpar0_no_drive_01_upwind,
                        sim_st_b001_fbpar0_no_drive_02_upwind,
                        sim_gs2_b001_fbpar0,
                           ],
                          [
                           "stella, zero drive",
                           "stella, zero drive, z & vpa upwind=0",
                           "stella, zero drive, z & vpa upwind=0.1",
                           "stella, zero drive, z & vpa upwind=0.2",
                           "GS2",
                           ],
                          IMAGE_DIR + "fbpar0_beta1e-3_zerodrive_scan_upwind",
                          sim_types=[
                                     "stella",
                                     "stella",
                                     "stella",
                                     "stella",
                                     "gs2",
                                     ],
                           plot_apar=True,
                           plot_bpar=False,
                           plot_format=".png"
                           )

def compare_omega_for_fbpar0_changing_streaming_and_drive():
    """ """
    make_comparison_plots([
                            sim_st_b001_fbpar0,
                            sim_st_b001_fbpar0_no_drive,
                            sim_st_b001_fbpar0_no_stream_no_mirror,
                            sim_st_b001_fbpar0_no_drive_no_stream_no_mirror,
                            sim_st_b001_fbpar0_no_drive_no_drifts
                         ],
                        [
                         "stella",
                         "stella, zero drive",
                         "stella, zero streaming & mirror",
                         "stella, zero drive, zero streaming & mirror",
                         "stella, zero drive, zero drifts",
                         ],
                        IMAGE_DIR + "fbpar0_beta1e-3_zerodrive_change_str_mirr",
                        sim_types=[
                                   "stella",
                                   "stella",
                                   "stella",
                                   "stella",
                                   "stella",
                                   ],
                         plot_apar=True,
                         plot_bpar=False,
                         plot_format=".png", show_fig=True
                         )

def compare_beta03():
    """ """
    make_comparison_plots([
                            sim_st_b03,
                            sim_gs2_b03,
                            sim_gs2_b03_hvr,
                         ],
                        [
                         "stella",
                         "gs2",
                         "gs2, higher vres",
                         ],
                        IMAGE_DIR + "beta=0.03",
                        sim_types=[
                                   "stella",
                                   "gs2",
                                   "gs2",
                                   ],
                         plot_apar=True,
                         plot_bpar=True,
                         plot_format=".png", show_fig=True
                         )

def compare_beta03_detailed():
    """ """

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)

    sim_gs2_b03_outnc = sim_gs2_b03 + ".out.nc"
    theta, gds2, gds21, gds22, bmag, gradpar = extract_data_from_ncdf(sim_gs2_b03_outnc,
                                    'theta', 'gds2', 'gds21', 'gds22', 'bmag', 'gradpar')
    ax1.plot(theta/np.pi, gds2, c="black")
    ax1.plot(theta/np.pi, gds21, ls="--", c="black")
    ax1.plot(theta/np.pi, gds22, ls="-.", c="black")
    ax2.plot(theta/np.pi, bmag, c="black")
    ax2.plot(theta/np.pi, gradpar, ls="--", c="black")

    sim_st_b03_outnc = sim_st_b03 + ".out.nc"
    theta, gds2, gds21, gds22, bmag, gradpar = extract_data_from_ncdf(sim_st_b03_outnc,
                                    'zed', 'gds2', 'gds21', 'gds22', 'bmag', 'gradpar')
    ax1.plot(theta/np.pi, gds2, c="red")
    ax1.plot(theta/np.pi, gds21, ls="--", c="red")
    ax1.plot(theta/np.pi, gds22, ls="-.", c="red")
    ax2.plot(theta/np.pi, bmag, c="red")
    ax2.plot(theta/np.pi, gradpar, ls="--", c="red")

    plt.show()

    make_comparison_plots([
                            sim_st_b03,
                            sim_gs2_b03,
                         ],
                        [
                         "stella",
                         "gs2",
                         ],
                        IMAGE_DIR + "beta=0.03",
                        sim_types=[
                                   "stella",
                                   "gs2",
                                   ],
                         plot_apar=True,
                         plot_bpar=True,
                         plot_format=".png", show_fig=False
                         )
    make_comparison_plots([
                            sim_st_b03_fbpar0,
                            sim_gs2_b03_fbpar0,
                         ],
                        [
                         "stella",
                         "gs2",
                         ],
                        IMAGE_DIR + "beta=0.03_fbpar0",
                        sim_types=[
                                   "stella",
                                   "gs2",
                                   ],
                         plot_apar=True,
                         plot_bpar=False,
                         plot_format=".png", show_fig=False
                         )
    make_comparison_plots([
                            sim_st_b03_fapar0,
                            sim_gs2_b03_fapar0,
                         ],
                        [
                         "stella",
                         "gs2",
                         ],
                        IMAGE_DIR + "beta=0.03_fapar0",
                        sim_types=[
                                   "stella",
                                   "gs2",
                                   ],
                         plot_apar=False,
                         plot_bpar=True,
                         plot_format=".png", show_fig=False
                         )
    return

def plot_g_for_fbpar0_different_terms_off():
    """Take a look at the distribrution function for
    fbpar=0 sims with beta=1e-3 and various terms on/off """

    stella_outnc_longname = sim_st_b001_fbpar0 + ".out.nc"
    stella_outnc_longname_no_mirror = sim_st_b001_fbpar0_no_mirror + ".out.nc"

    plot_gmvus(stella_outnc_longname)
    plot_gzvs(stella_outnc_longname)
    print("Now no mirror:")
    plot_gmvus(stella_outnc_longname_no_mirror)
    plot_gzvs(stella_outnc_longname_no_mirror)


    return

def plot_gzvs_for_fbpar0():
    """Take a look at the distribrution function for
    fbpar=0 sims """

    stella_sim_longname = "stella_beta0.010_fbpar0/input"
    stella_outnc_longname = stella_sim_longname + ".out.nc"
    plot_gzvs(stella_outnc_longname)


    return

def make_comparison_plots_many(stella_sim_longnames, gs2_sim_longnames,
                                beta_strs, prefix, plot_apar = False,
                                plot_bpar = False, plot_format=".png"):
    """ """

    for beta_idx in range(0, len(beta_strs)):
        stella_sim_longname = stella_sim_longnames[beta_idx]
        gs2_sim_longname = gs2_sim_longnames[beta_idx]
        beta_str = beta_strs[beta_idx]
        make_comparison_plots([
                               stella_sim_longname,
                               gs2_sim_longname,
                               ],
                              [
                               "stella",
                               "GS2",
                               ],
                              IMAGE_DIR + prefix + "beta_" + beta_str,
                              sim_types=[
                                         "stella",
                                         "gs2",
                                         ],
                               plot_apar=plot_apar,
                               plot_bpar=plot_bpar,
                               plot_format=plot_format
                               )

def plot_fapar_fbpar_on():
    """ """
    ## Beta scan
    stella_sim_longnames = [
                            sim_st_b0,
                            sim_st_b005,
                            sim_st_b01,
                            sim_st_b015,
                            sim_st_b02,
                            sim_st_b025,
                            sim_st_b03
                            ]
    stella_beta_vals = [
                        0.,
                        0.005,
                        0.01,
                        0.015,
                        0.02,
                        0.025,
                        0.03,
                        ]
    stella_labels = [
                        "beta=0.",
                        "beta=0.005",
                        "beta=0.01",
                        "beta=0.015",
                        "beta=0.02",
                        "beta=0.025",
                        "beta=0.03",
                        ]
    make_comparison_plots(stella_sim_longnames,
                          stella_labels,
                          "images/omega_beta_scan/fapar1_fbpar1_beta_scan",
                          plot_apar=True, plot_bpar=True, plot_format=".eps")

    make_beta_scan_plots(stella_sim_longnames,
                            [],
                         stella_beta_vals,
                         IMAGE_DIR + "test_cbc_beta_scan",
                         gs2_pickle=pickle_gs2)

    return

def plot_fapar0():
    """ """
    stella_sim_longnames = [
                            sim_st_b00001_fapar0,
                            sim_st_b001_fapar0,
                            sim_st_b002_fapar0,
                            sim_st_b01_fapar0,
                            ]
    gs2_sim_longnames = [
                         sim_gs2_b00001_fapar0,
                         sim_gs2_b001_fapar0,
                         sim_gs2_b002_fapar0,
                         sim_gs2_b01_fapar0
                        ]
    beta_strs = [
                 "0.00001",
                 "0.001",
                 "0.002",
                 "0.01"
                 ]
    # make_comparison_plots_many(stella_sim_longnames,
    #                            gs2_sim_longnames,
    #                            beta_strs, "fapar=0/", plot_apar=False, plot_bpar=True)

    make_beta_scan_plots(stella_sim_longnames,
                         gs2_sim_longnames,
                         beta_strs,
                         IMAGE_DIR + "test_cbc_beta_scan_fapar0",
                         )

    return

def plot_fbpar0():
    """ """
    stella_sim_longnames = [
                            sim_st_b00001_fbpar0,
                            sim_st_b00005_fbpar0,
                            sim_st_b0001_fbpar0,
                            sim_st_b0003_fbpar0,
                            sim_st_b0006_fbpar0,
                            #sim_st_b001_fbpar0,
                            sim_st_b0015_fbpar0,
                            sim_st_b002_fbpar0,
                            sim_st_b003_fbpar0,
                            sim_st_b004_fbpar0,
                            sim_st_b005_fbpar0,
                            sim_st_b01_fbpar0
                            ]

    gs2_sim_longnames = [
                            sim_gs2_b00001_fbpar0,
                            sim_gs2_b001_fbpar0,
                            sim_gs2_b001_fbpar0,
                            sim_gs2_b001_fbpar0,
                            sim_gs2_b001_fbpar0,
                            #sim_gs2_b001_fbpar0,
                            sim_gs2_b001_fbpar0,
                            sim_gs2_b002_fbpar0,
                            sim_gs2_b003_fbpar0,
                            sim_gs2_b004_fbpar0,
                            sim_gs2_b005_fbpar0,
                            sim_gs2_b01_fbpar0
                        ]
    beta_strs = [
                 "0.00001",
                 "0.00005",
                 "0.0001",
                 "0.0003",
                 "0.0006",
                 "0.001",
                 "0.0015",
                 "0.002",
                 "0.003",
                 "0.004",
                 "0.005",
                 "0.01"
                 ]
    # make_comparison_plots_many(stella_sim_longnames,
    #                            gs2_sim_longnames,
    #                            beta_strs, "fbpar=0/", plot_apar=True, plot_bpar=False)
    make_beta_scan_plots(stella_sim_longnames,
                         gs2_sim_longnames,
                         beta_strs,
                         IMAGE_DIR + "test_cbc_beta_scan_fbpar0",
                             )
    return

def plot_fbpar0_with_gs2_pickle():
    """ """
    stella_sim_longnames = [
                            sim_st_b00001_fbpar0,
                            sim_st_b00005_fbpar0,
                            sim_st_b0001_fbpar0,
                            sim_st_b0003_fbpar0,
                            sim_st_b0006_fbpar0,
                            #sim_st_b001_fbpar0,
                            sim_st_b0015_fbpar0,
                            sim_st_b002_fbpar0,
                            sim_st_b003_fbpar0,
                            sim_st_b004_fbpar0,
                            sim_st_b005_fbpar0,
                            sim_st_b01_fbpar0
                            ]

    gs2_sim_longnames = [
                        ]
    beta_strs = [
                 "0.00001",
                 "0.00005",
                 "0.0001",
                 "0.0003",
                 "0.0006",
                 "0.001",
                 "0.0015",
                 "0.002",
                 "0.003",
                 "0.004",
                 "0.005",
                 "0.01"
                 ]
    # make_comparison_plots_many(stella_sim_longnames,
    #                            gs2_sim_longnames,
    #                            beta_strs, "fbpar=0/", plot_apar=True, plot_bpar=False)
    make_beta_scan_plots(stella_sim_longnames,
                         [],
                         beta_strs,
                         IMAGE_DIR + "test_cbc_beta_scan_fbpar0",
                         gs2_pickle=pickle_gs2_fbpar0     )
    return

def compare_get_fields_subroutines():
    """ """
    make_comparison_plots(
                    [
                    sim_st_gfvmulo_nstep1000,
                    sim_st_gf_nstep1000,
                    ],
                    [
                    "get_fields_vmulo",
                    "get_fields",
                    ],
                    IMAGE_DIR ,
                    sim_types = [
                            "stella",
                            "stella",
                                ],
                     plot_apar=True,
                     plot_bpar=True,
                     plot_format=".png", show_fig=True
                    )


    return

def make_all_plots():
    """ """
    plot_fbpar0()
    plot_fapar_fbpar_on()
    plot_fapar0()
    compare_omega_for_fbpar0_zero_drive()
    compare_omega_for_fbpar0_zero_drive_change_upwind()
    analyse_fbpar0_beta0001_results()

    return

if __name__ == "__main__":
    ## Compare

    #analyse_fbpar0_results()
    # plot_beta_scans()
    #analyse_results_for_poster()
    # plot_geometry()
    # make_low_beta_fbpar0_plots()
    # analyse_fapar0_results()
    # plot_gvmus_for_fbpar0()
    # analyse_fapar0_changing_vpares()
    #make_all_plots()
    #plot_gzvs_for_fbpar0()
    # plot_fapar0()
    # plot_fbpar0_with_gs2_pickle()
    # plot_fapar_fbpar_on()
    #plot_geometry()
    #compare_beta03()
    #compare_beta03_detailed()
    # compare_omega_for_fbpar0_changing_streaming_and_drive()
    #analyse_fbpar0_beta0001_results()
    #plot_fbpar0_beta0001_equal_masses()
    #plot_g_for_fbpar0_different_terms_off()
    #compare_omega_for_fbpar0_different_terms_off()
    # compare_omega_for_fbpar0_zero_drive()
    # compare_omega_for_fbpar0_zero_drive_change_upwind()
    #compare_omega_for_fbpar0_changing_streaming_and_drive()
    compare_get_fields_subroutines()
