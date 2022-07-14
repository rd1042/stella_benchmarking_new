"""See what happens to KBM in CBC at beta=4% in stella/GS2 as we vary:
 - dt
 - space-centering
 - time-centering
 - vpa-centering
 - v-space extent/resolution
 - theta extent/resolution """

import sys
import pickle
import numpy as np
import glob
import re
import matplotlib.pyplot as plt
from itertools import cycle
import xarray as xr

sys.path.append("../postprocessing_tools")
import helper_linear_sims as help_lin
from helper_ncdf import view_ncdf_variables
from extract_sim_data import get_omega_data, get_phiz_data, get_aparz_data, get_bparz_data
from plotting_helper import plot_phi_z_for_sim, plot_apar_z_for_sim, plot_bpar_z_for_sim
from helper_ncdf import view_ncdf_variables
import plot_2d_utils as plot2dutils

comparison_folder = "sims/beta_0.04000_investigation/"

def make_comparison(outnc_longnames, gamma_buffer = 0.01, omega_buffer = 0.01,lw=3):
    """ """
    sim_shortnames = []
    tvals_list = []
    tend_list = []
    zvals_list = []
    gamma_list = []
    omega_list = []
    gamma_final_list = []
    omega_final_list = []
    phiz_list = []
    aparz_list = []
    bparz_list = []
    for outnc_longname in outnc_longnames:
        ## For each sim, we want to compare:
        # Mode structures
        # Omega(t) (are we converged?)
        # Omega (print to screen?)
        try:
            sim_shortname = re.split("/", outnc_longname)[-1]
            sim_shortname = re.split(".out.nc", sim_shortname)[0]
            sim_longname = re.split(".out.nc", outnc_longname)[0]
            # view_ncdf_variables(outnc_longname)
            if "gs2" in sim_shortname:
                sim_type = "gs2"
                # print("sim_name: ", sim_shortname)
                # view_ncdf_variables(outnc_longname)
            else:
                sim_type = "stella"
                # print("sim_name: ", sim_shortname)
                # data = xr.open_dataset(outnc_longname)
                # print("data.keys = ", data.keys())
                # sys.exit()
            time, freqom_final, gammaom_final, freqom, gammaom, gamma_stable = get_omega_data(sim_longname, sim_type)
            #try:
            z, real_phi, imag_phi = get_phiz_data(sim_longname, sim_type)
            abs_phi = help_lin.get_abs(real_phi, imag_phi)
            try:
                z, real_apar, imag_apar = get_aparz_data(sim_longname, sim_type)
                abs_apar = help_lin.get_abs(real_apar, imag_apar)
            except None:
                print("None")
            try:
                z, real_bpar, imag_bpar = get_bparz_data(sim_longname, sim_type)
                abs_bpar = help_lin.get_abs(real_bpar, imag_bpar)
            except None:
                print("None")
            max_abs_phi = np.max(abs_phi)
            sim_shortnames.append(sim_shortname)
            tvals_list.append(time)
            tend_list.append(time[-1])
            zvals_list.append(z)
            gamma_list.append(gamma_stable)
            omega_list.append(freqom)
            gamma_final_list.append(gamma_stable[-1])
            omega_final_list.append(freqom[-1])
            phiz_list.append(abs_phi/max_abs_phi)
            aparz_list.append(abs_apar/max_abs_phi)
            bparz_list.append(abs_bpar/max_abs_phi)
            print(sim_shortname, freqom[-1], gamma_stable[-1])
        except (IndexError, FileNotFoundError):
            # Probably caused by an incomplete simulation.
            pass

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    omega_linestyles=cycle(["-","--","-."])
    gamma_linestyles=cycle(["-","--","-."])
    for sim_idx, sim_shortname in enumerate(sim_shortnames):
        ax1.plot(tvals_list[sim_idx]-tend_list[sim_idx], omega_list[sim_idx], label=sim_shortname, ls=next(omega_linestyles))
        ax2.plot(tvals_list[sim_idx]-tend_list[sim_idx], gamma_list[sim_idx], label=sim_shortname, ls=next(gamma_linestyles))

    max_omega_final = np.max(np.array(omega_final_list))
    min_omega_final = np.min(np.array(omega_final_list))
    max_gamma_final = np.max(np.array(gamma_final_list))
    min_gamma_final = np.min(np.array(gamma_final_list))

    ax1.legend(loc="best")
    ax2.set_xlabel(r"$t$")
    ax1.set_ylabel(r"$\omega$")
    ax2.set_ylabel(r"$\gamma$")
    ax1.set_ylim(min_omega_final-omega_buffer, max_omega_final+omega_buffer)
    ax2.set_ylim(min_gamma_final-gamma_buffer, max_gamma_final+gamma_buffer)

    for ax in [ax1, ax2]:
        ax.grid(True)

    plt.show()

    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312, sharex=ax1)
    ax3 = fig.add_subplot(313, sharex=ax1)
    phi_linestyles=cycle(["-","--","-."])
    apar_linestyles=cycle(["-","--","-."])
    bpar_linestyles=cycle(["-","--","-."])
    for sim_idx, sim_shortname in enumerate(sim_shortnames):
        ax1.plot(zvals_list[sim_idx]/np.pi, phiz_list[sim_idx], label=sim_shortname, ls=next(phi_linestyles))
        ax2.plot(zvals_list[sim_idx]/np.pi, aparz_list[sim_idx], label=sim_shortname, ls=next(apar_linestyles))
        ax3.plot(zvals_list[sim_idx]/np.pi, bparz_list[sim_idx], label=sim_shortname, ls=next(bpar_linestyles))

    ax1.legend(loc="best")
    ax3.set_xlabel(r"$z/\pi$")
    ax1.set_ylabel(r"$\phi$")
    ax2.set_ylabel(r"$A_\parallel$")
    ax3.set_ylabel(r"$B_\parallel$")
    # ax1.set_ylim(min_omega_final-omega_buffer, max_omega_final+omega_buffer)
    # ax2.set_ylim(min_gamma_final-gamma_buffer, max_gamma_final+gamma_buffer)

    for ax in [ax1, ax2, ax3]:
        ax.grid(True)

    plt.show()
    return

def compare_all_sims():
    """ """
    # See all files except blowup
    outnc_longnames = glob.glob(comparison_folder + "[!blowup]*/*.out.nc")

    make_comparison(outnc_longnames)
    return

def compare_ntheta32_sims():
    """ """
    # See GS2 and stella ntheta=32
    outnc_longnames_1 = glob.glob(comparison_folder + "gs2/*.out.nc")
    outnc_longnames_2 = glob.glob(comparison_folder + "stella_str_impl_ntheta32_mirror_explicit/*.out.nc")
    outnc_lists = [outnc_longnames_1, outnc_longnames_2]
    # Flatten into a 1D list
    outnc_longnames =  [element for sublist in outnc_lists for element in sublist]
    make_comparison(outnc_longnames)
    return

def compare_ntheta64_sims():
    """ """
    # See GS2 and stella ntheta=32
    outnc_longnames_1 = glob.glob(comparison_folder + "gs2/*.out.nc")
    outnc_longnames_2 = glob.glob(comparison_folder + "stella_str_impl_ntheta64/*.out.nc")
    outnc_lists = [outnc_longnames_1, outnc_longnames_2]
    # Flatten into a 1D list
    outnc_longnames =  [element for sublist in outnc_lists for element in sublist]
    make_comparison(outnc_longnames)
    return

def compare_ntheta32_vs_64_sims():
    """ """
    # See GS2 and stella ntheta=32 and ntheta=64 (fiducial points only)
    outnc_longnames_1 = glob.glob(comparison_folder + "gs2/*.out.nc")
    outnc_longnames_2 = glob.glob(comparison_folder + "stella_str_impl_ntheta32_mirror_explicit/stella_str_impl_ntheta32.out.nc")
    outnc_longnames_3 = glob.glob(comparison_folder + "stella_str_impl_ntheta64/stella_str_impl.out.nc")
    outnc_lists = [outnc_longnames_1, outnc_longnames_2, outnc_longnames_3]
    # Flatten into a 1D list
    outnc_longnames =  [element for sublist in outnc_lists for element in sublist]
    make_comparison(outnc_longnames)
    return

def compare_blowup():
    """ """
    outnc_longnames = glob.glob(comparison_folder + "blowup/*.out.nc")
    make_comparison(outnc_longnames, gamma_buffer = 3, omega_buffer = 80)
    return

def compate_dibrn_fn_for_2_sims(outnc_longnames, sim_labels):
    """ """

    def make_gvmus_plots():
        """ """
        fig = plt.figure(figsize=[12, 12])
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)
        for mu_idx, mu_val in enumerate(mu_0):
            ax1.plot(vpa_0, gvmus_tfinal_ion_0[mu_idx,:], label="mu={:.2f}".format(float(mu_val)))
            ax2.plot(vpa_0, gvmus_tfinal_electron_0[mu_idx,:], label="mu={:.2f}".format(float(mu_val)))
        for mu_idx, mu_val in enumerate(mu_1):
            ax3.plot(vpa_1, gvmus_tfinal_ion_1[mu_idx,:], label="mu={:.2f}".format(float(mu_val)))
            ax4.plot(vpa_1, gvmus_tfinal_electron_1[mu_idx,:], label="mu={:.2f}".format(float(mu_val)))

        for ax in [ax1, ax2, ax3, ax4]:
            ax.grid(True)
            ax.set_xlabel(r"$v_\parallel$")
            ax.legend(loc="best")
        ax1.set_title(sim_labels[0])
        ax2.set_title(sim_labels[0])
        ax3.set_title(sim_labels[1])
        ax4.set_title(sim_labels[1])
        ax1.set_ylabel("gvmus(ion)")
        ax2.set_ylabel("gvmus(electron)")
        ax3.set_ylabel("gvmus(ion)")
        ax4.set_ylabel("gvmus(electron)")
        plt.show()

        fig = plt.figure(figsize=[12, 12])
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)
        for vpa_idx, vpa_val in enumerate(vpa_0):
            ax1.plot(mu_0, gvmus_tfinal_ion_0[:,vpa_idx], label="vpa={:.2f}".format(float(vpa_val)))
            ax2.plot(mu_0, gvmus_tfinal_electron_0[:,vpa_idx], label="vpa={:.2f}".format(float(vpa_val)))
        for vpa_idx, vpa_val in enumerate(vpa_1):
            ax3.plot(mu_1, gvmus_tfinal_ion_1[:,vpa_idx], label="vpa={:.2f}".format(float(vpa_val)))
            ax4.plot(mu_1, gvmus_tfinal_electron_1[:,vpa_idx], label="vpa={:.2f}".format(float(vpa_val)))

        for ax in [ax1, ax2, ax3, ax4]:
            ax.grid(True)
            ax.set_xlabel(r"$\mu$")
            ax.legend(loc="best")
        ax1.set_title(sim_labels[0])
        ax2.set_title(sim_labels[0])
        ax3.set_title(sim_labels[1])
        ax4.set_title(sim_labels[1])
        ax1.set_ylabel("gvmus(ion)")
        ax2.set_ylabel("gvmus(electron)")
        ax3.set_ylabel("gvmus(ion)")
        ax4.set_ylabel("gvmus(electron)")
        plt.show()


        mu_meshgrid, vpa_meshgrid, [gion_0_meshgrid, gelectron_0_meshgrid] = plot2dutils.uniquearrays2meshgrids(mu_0,
                    vpa_0, [gvmus_tfinal_ion_0, gvmus_tfinal_electron_0], n_xpoints=400, n_ypoints=400)

        fig = plt.figure(figsize=[12, 12])
        plot_dim_x = 0.36
        plot_dim_y = plot_dim_x
        col1_lhs = 0.05
        col2_lhs = 0.53
        col1_cb_lhs = 0.42
        col2_cb_lhs = 0.93
        row1_bottom = 0.54
        row2_bottom = 0.08
        cb_width = 0.02
        ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
        ax2 = fig.add_axes([col2_lhs, row1_bottom, plot_dim_x, plot_dim_y])
        ax3 = fig.add_axes([col1_lhs, row2_bottom, plot_dim_x, plot_dim_y])
        ax4 = fig.add_axes([col2_lhs, row2_bottom, plot_dim_x, plot_dim_y])
        cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
        cbax2 = fig.add_axes([col2_cb_lhs, row1_bottom, cb_width, plot_dim_y])
        cbax3 = fig.add_axes([col1_cb_lhs, row2_bottom, cb_width, plot_dim_y])
        cbax4 = fig.add_axes([col2_cb_lhs, row2_bottom, cb_width, plot_dim_y])

        g_levels = 40
        # Filled contour plot for growth rate
        gion_0_contours = ax1.contourf(mu_meshgrid, vpa_meshgrid, gion_0_meshgrid,
                                          g_levels, cmap="inferno", extend="max")
        fig.colorbar(gion_0_contours, cax=cbax1)
        gelectron_0_contours = ax2.contourf(mu_meshgrid, vpa_meshgrid, gelectron_0_meshgrid,
                                          g_levels, cmap="inferno", extend="max")
        fig.colorbar(gelectron_0_contours, cax=cbax2)

        mu_meshgrid, vpa_meshgrid, [gion_1_meshgrid, gelectron_1_meshgrid] = plot2dutils.uniquearrays2meshgrids(mu_1,
                    vpa_1, [gvmus_tfinal_ion_1, gvmus_tfinal_electron_1], n_xpoints=400, n_ypoints=400)

        gion_1_contours = ax3.contourf(mu_meshgrid, vpa_meshgrid, gion_1_meshgrid,
                                          g_levels, cmap="inferno", extend="max")
        fig.colorbar(gion_1_contours, cax=cbax3)
        gelectron_1_contours = ax4.contourf(mu_meshgrid, vpa_meshgrid, gelectron_1_meshgrid,
                                          g_levels, cmap="inferno", extend="max")
        fig.colorbar(gelectron_1_contours, cax=cbax4)

        for ax in [ax1, ax2, ax3, ax4]:
            ax.set_xlabel(r"$\mu$")
            ax.set_ylabel(r"$v_\parallel$")

        ax1.set_title(sim_labels[0] + ", ion")
        ax2.set_title(sim_labels[0] + ", electron")
        ax3.set_title(sim_labels[1] + ", ion")
        ax4.set_title(sim_labels[1] + ", electron")

        plt.show()

        return

    def make_gvzs_plots():
        """ """
        fig = plt.figure(figsize=[12, 12])
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)
        for z_idx, z_val in enumerate(z_0):
            ax1.plot(vpa_0, gzvs_tfinal_ion_0[:,z_idx], label="z={:.2f}".format(float(z_val)))
            ax2.plot(vpa_0, gzvs_tfinal_electron_0[:,z_idx], label="z={:.2f}".format(float(z_val)))
        for z_idx, z_val in enumerate(z_1):
            ax3.plot(vpa_1, gzvs_tfinal_ion_1[:,z_idx], label="z={:.2f}".format(float(z_val)))
            ax4.plot(vpa_1, gzvs_tfinal_electron_1[:,z_idx], label="z={:.2f}".format(float(z_val)))

        for ax in [ax1, ax2, ax3, ax4]:
            ax.grid(True)
            ax.set_xlabel(r"$v_\parallel$")
            ax.legend(loc="best")
        ax1.set_title(sim_labels[0])
        ax2.set_title(sim_labels[0])
        ax3.set_title(sim_labels[1])
        ax4.set_title(sim_labels[1])
        ax1.set_ylabel("gzvs(ion)")
        ax2.set_ylabel("gzvs(electron)")
        ax3.set_ylabel("gzvs(ion)")
        ax4.set_ylabel("gzvs(electron)")
        plt.show()

        fig = plt.figure(figsize=[12, 12])
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)
        for vpa_idx, vpa_val in enumerate(vpa_0):
            ax1.plot(z_0/np.pi, gzvs_tfinal_ion_0[vpa_idx,:], label="vpa={:.2f}".format(float(vpa_val)))
            ax2.plot(z_0/np.pi, gzvs_tfinal_electron_0[vpa_idx,:], label="vpa={:.2f}".format(float(vpa_val)))
        for vpa_idx, vpa_val in enumerate(vpa_1):
            ax3.plot(z_1/np.pi, gzvs_tfinal_ion_1[vpa_idx,:], label="vpa={:.2f}".format(float(vpa_val)))
            ax4.plot(z_1/np.pi, gzvs_tfinal_electron_1[vpa_idx,:], label="vpa={:.2f}".format(float(vpa_val)))

        for ax in [ax1, ax2, ax3, ax4]:
            ax.grid(True)
            ax.set_xlabel(r"$z/\pi$")
            ax.legend(loc="best")
        ax1.set_title(sim_labels[0])
        ax2.set_title(sim_labels[0])
        ax3.set_title(sim_labels[1])
        ax4.set_title(sim_labels[1])
        ax1.set_ylabel("gzvs(ion)")
        ax2.set_ylabel("gzvs(electron)")
        ax3.set_ylabel("gzvs(ion)")
        ax4.set_ylabel("gzvs(electron)")
        plt.show()

        z_meshgrid, vpa_meshgrid, [gion_0_meshgrid, gelectron_0_meshgrid] = plot2dutils.uniquearrays2meshgrids(z_0,
                    vpa_0, [gzvs_tfinal_ion_0, gzvs_tfinal_electron_0], n_xpoints=400, n_ypoints=400)

        fig = plt.figure(figsize=[12, 12])
        plot_dim_x = 0.36
        plot_dim_y = plot_dim_x
        col1_lhs = 0.05
        col2_lhs = 0.53
        col1_cb_lhs = 0.42
        col2_cb_lhs = 0.93
        row1_bottom = 0.54
        row2_bottom = 0.08
        cb_width = 0.02
        ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
        ax2 = fig.add_axes([col2_lhs, row1_bottom, plot_dim_x, plot_dim_y])
        ax3 = fig.add_axes([col1_lhs, row2_bottom, plot_dim_x, plot_dim_y])
        ax4 = fig.add_axes([col2_lhs, row2_bottom, plot_dim_x, plot_dim_y])
        cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])
        cbax2 = fig.add_axes([col2_cb_lhs, row1_bottom, cb_width, plot_dim_y])
        cbax3 = fig.add_axes([col1_cb_lhs, row2_bottom, cb_width, plot_dim_y])
        cbax4 = fig.add_axes([col2_cb_lhs, row2_bottom, cb_width, plot_dim_y])

        g_levels = 40
        # Filled contour plot for growth rate
        gion_0_contours = ax1.contourf(z_meshgrid, vpa_meshgrid, gion_0_meshgrid,
                                          g_levels, cmap="inferno", extend="max")
        fig.colorbar(gion_0_contours, cax=cbax1)
        gelectron_0_contours = ax2.contourf(z_meshgrid, vpa_meshgrid, gelectron_0_meshgrid,
                                          g_levels, cmap="inferno", extend="max")
        fig.colorbar(gelectron_0_contours, cax=cbax2)

        z_meshgrid, vpa_meshgrid, [gion_1_meshgrid, gelectron_1_meshgrid] = plot2dutils.uniquearrays2meshgrids(z_1,
                    vpa_1, [gzvs_tfinal_ion_1, gzvs_tfinal_electron_1], n_xpoints=400, n_ypoints=400)

        gion_1_contours = ax3.contourf(z_meshgrid, vpa_meshgrid, gion_1_meshgrid,
                                          g_levels, cmap="inferno", extend="max")
        fig.colorbar(gion_1_contours, cax=cbax3)
        gelectron_1_contours = ax4.contourf(z_meshgrid, vpa_meshgrid, gelectron_1_meshgrid,
                                          g_levels, cmap="inferno", extend="max")
        fig.colorbar(gelectron_1_contours, cax=cbax4)

        for ax in [ax1, ax2, ax3, ax4]:
            ax.set_xlabel(r"$z$")
            ax.set_ylabel(r"$v_\parallel$")

        ax1.set_title(sim_labels[0] + ", ion")
        ax2.set_title(sim_labels[0] + ", electron")
        ax3.set_title(sim_labels[1] + ", ion")
        ax4.set_title(sim_labels[1] + ", electron")
        plt.show()
        return

    data_0 = xr.open_dataset(outnc_longnames[0])
    data_1 = xr.open_dataset(outnc_longnames[1])
    # print("data.keys = ", data.keys())
    t_1 = data_1.t
    z_1 = data_1.zed
    vpa_1 = data_1.vpa
    mu_1 = data_1.mu
    gzvs_1 = data_1.gzvs
    gvmus_1 = data_1.gvmus
    t_0 = data_0.t
    z_0 = data_0.zed
    vpa_0 = data_0.vpa
    mu_0 = data_0.mu
    gzvs_0 = data_0.gzvs
    gvmus_0 = data_0.gvmus

    gzvs_tfinal_ion_1 = gzvs_1[-1,0,:,:,0]
    gzvs_tfinal_electron_1 = gzvs_1[-1,1,:,:,0]
    gvmus_tfinal_ion_1 = gvmus_1[-1,0,:,:]
    gvmus_tfinal_electron_1 = gvmus_1[-1,1,:,:]
    gzvs_tfinal_ion_0 = gzvs_0[-1,0,:,:,0]
    gzvs_tfinal_electron_0 = gzvs_0[-1,1,:,:,0]
    gvmus_tfinal_ion_0 = gvmus_0[-1,0,:,:]
    gvmus_tfinal_electron_0 = gvmus_0[-1,1,:,:]

    make_gvmus_plots()
    make_gvzs_plots()
    #print("max(abs(gvmus_tfinal_ion - gvmus_tfinal_electron)) = ", np.max(abs(gvmus_tfinal_ion - gvmus_tfinal_electron)))

    return

def examine_moments_for_blowup():
    """ """
    outnc_longname_fapar0_fbpar0 = "sims/beta_0.04000_investigation/blowup/stella_str_impl_mirror_impl_dt5E-3_longer_run_time_and_moments_fapar0_fbpar0.out.nc"
    outnc_longname_fapar0_fbpar1 = "sims/beta_0.04000_investigation/blowup/stella_str_impl_mirror_impl_dt5E-3_longer_run_time_and_moments_fapar0.out.nc"
    outnc_longname_fapar1_fbpar0 = "sims/beta_0.04000_investigation/blowup/stella_str_impl_mirror_impl_dt5E-3_longer_run_time_and_moments_fbpar0.out.nc"
    outnc_longname_fapar1_fbpar1_stable = "sims/beta_0.04000_investigation/stella_str_impl_ntheta32_mirror_implicit/stella_ntheta32_str_impl_mirror_impl_dt5E-3_longer_run_time_and_moments.out.nc"
    outnc_longname_fapar1_fbpar1_unstable = "sims/beta_0.04000_investigation/blowup/stella_str_impl_mirror_impl_dt5E-3_longer_run_time_and_moments.out.nc"

    # compate_dibrn_fn_for_2_sims([outnc_longname_fapar1_fbpar1_stable, outnc_longname_fapar1_fbpar1_unstable],
    #                             [r"$n_z=32$(stable)", r"$n_z=64$(unstable)"])
    compate_dibrn_fn_for_2_sims([outnc_longname_fapar0_fbpar0, outnc_longname_fapar1_fbpar1_unstable],
                                [r"fapar=fbpar=0 (stable)", r"$n_z=64$ (unstable)"])
    # compate_dibrn_fn_for_2_sims([outnc_longname_fapar1_fbpar0, outnc_longname_fapar1_fbpar1_unstable],
    #                             [r"fapar=1,fbpar=0 (stable)", r"$n_z=64$ (unstable)"])
    # compate_dibrn_fn_for_2_sims([outnc_longname_fapar0_fbpar1, outnc_longname_fapar1_fbpar1_unstable],
    #                             [r"fapar=0,fbpar=1 (unstable)", r"$n_z=64$ (unstable)"])

    return

if __name__ == "__main__":
    print("Hello world")
    # compare_blowup()
    # compare_ntheta32_vs_64_sims()
    # compare_ntheta32_sims()
    # compare_ntheta64_sims()
    # compare_all_sims()
    examine_moments_for_blowup()
