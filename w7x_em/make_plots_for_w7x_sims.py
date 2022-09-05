""" """

import sys
sys.path.append("../postprocessing_tools")
from helper_ncdf_new import view_ncdf_variables, extract_data_from_ncdf
import matplotlib.pyplot as plt
import numpy as np
import helper_linear_sims as help_lin
from extract_sim_data import get_omega_data, get_phiz_data, get_aparz_data, get_bparz_data
from plotting_helper import plot_phi_z_for_sim, plot_apar_z_for_sim, plot_bpar_z_for_sim
import re

def get_gradient_from_outnc_longname(outnc_longname):
    """ """

    sim_shortname = re.split("/", outnc_longname)[-1]
    sim_shortname = re.split(".out.nc", sim_shortname)[0]
    gradients_str = re.split("gradients_", sim_shortname)[-1]
    gradients_str = re.split("_", gradients_str)[0]
    return float(gradients_str)


def make_omega_time_plot_for_stella(sim_longname, time, freqom, gammaom, gamma_stable):
    """ """
    fig = plt.figure(figsize=(12,8))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    ax1.plot(time, freqom)
    ax2.plot(time, gammaom, label="from .omega")
    ax2.plot(time, gamma_stable, label="from phi2")
    ax1.set_ylabel(r"$\omega$")
    ax2.set_ylabel(r"$\gamma$")
    ax2.set_xlabel(r"$t$")

    ## Set axis limits; convergence should be better than 0.1, so use this
    tolerance=0.1
    if np.isfinite(freqom[-1]):
        ax1.set_ylim(freqom[-1]-tolerance, freqom[-1]+tolerance)
    if np.isfinite(gamma_stable[-1]):
        ax2.set_ylim(gamma_stable[-1]-tolerance, gamma_stable[-1]+tolerance)

    for ax in [ax1, ax2]:
        ax.grid(True)

    save_name = sim_longname + "_omega_t.png"
    plt.savefig(save_name)
    plt.close()

    return

def make_field_z_plot_for_stella(sim_longname, z, abs_phi, abs_apar, abs_bpar):
    """ """
    fig = plt.figure(figsize=(12,8))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312, sharex=ax1)
    ax3 = fig.add_subplot(313, sharex=ax1)
    norm_val = np.max(abs_phi)
    ax1.plot(z/np.pi, abs_phi/norm_val)
    ax2.plot(z/np.pi, abs_apar/norm_val)
    ax3.plot(z/np.pi, abs_bpar/norm_val)
    ax1.set_ylabel(r"$\vert \phi \vert$")
    ax2.set_ylabel(r"$\vert A_\parallel \vert$")
    ax3.set_ylabel(r"$\vert B_\parallel \vert$")
    ax3.set_xlabel(r"$z/\pi$")

    for ax in [ax1, ax2, ax3]:
        ax.grid(True)
    plt.tight_layout()
    save_name = sim_longname + "_fields_z.png"
    plt.savefig(save_name)
    plt.close()
    return


def postprocess_folder_stella(outnc_longnames, param_scanned="gradients",
                            save_as_pickle=False, make_plot=False):
    """For a folder, get Omega and beta, and save to pickle"""

    gradient_longlist = []
    growth_rate_longlist = []
    freq_longlist = []
    z_longlist = []
    abs_phi_longlist = []
    abs_apar_longlist = []
    abs_bpar_longlist = []

    for outnc_longname in outnc_longnames:
        sim_shortname = re.split("/", outnc_longname)[-1]
        sim_shortname = re.split(".out.nc", sim_shortname)[0]
        sim_longname = re.split(".out.nc", outnc_longname)[0]
        # view_ncdf_variables(outnc_longname)
        ## Get beta
        gradient_longlist.append(get_gradient_from_outnc_longname(outnc_longname))
        time, freqom_final, gammaom_final, freqom, gammaom, gamma_stable = get_omega_data(sim_longname, "stella")
        #try:
        z, real_phi, imag_phi = get_phiz_data(sim_longname, "stella")
        abs_phi = help_lin.get_abs(real_phi, imag_phi)
        try:
            z, real_apar, imag_apar = get_aparz_data(sim_longname, "stella")
            abs_apar = help_lin.get_abs(real_apar, imag_apar)
        except None:
            print("None")
        try:
            z, real_bpar, imag_bpar = get_bparz_data(sim_longname, "stella")
            abs_bpar = help_lin.get_abs(real_bpar, imag_bpar)
        except None:
            print("None")

        if make_plot:
            make_omega_time_plot_for_stella(sim_longname, time, freqom, gammaom, gamma_stable)
            make_field_z_plot_for_stella(sim_longname, z, abs_phi, abs_apar, abs_bpar)
        z_longlist.append(z)
        abs_phi_longlist.append(abs_phi)
        abs_apar_longlist.append(abs_apar)
        abs_bpar_longlist.append(abs_bpar)
        # (fitted_growth_rate, growth_rate_error,
        #  converged_frequency, freq_error) = help_lin.calculate_omega_and_plot_for_single(outnc_longname,
        #                 save_path=(folder_longname + "/"),
        #                 plot_growth_rates=True,
        #                 figtitle=sim_shortname)

        # growth_rate_longlist.append(gamma_stable[-1])
        growth_rate_longlist.append(gammaom_final)
        freq_longlist.append(freqom_final)

    ## Sort beta and Omega
    gradient_vals = np.array(gradient_longlist)
    growth_rate_vals = np.array(growth_rate_longlist)
    freq_vals = np.array(freq_longlist)
    sort_idxs = np.argsort(gradient_vals)
    gradient_vals = gradient_vals[sort_idxs]
    growth_rate_vals = growth_rate_vals[sort_idxs]
    freq_vals = freq_vals[sort_idxs]

    ## Save to a pickle
    if save_as_pickle:
        pickle_name = folder_longname + "/beta_gamma_omega.pickle"
        myfile = open(pickle_name, "wb")
        pickle.dump([gradient_vals, growth_rate_vals, freq_vals], myfile)
        myfile.close()

        return
    else:
        return (gradient_longlist, z_longlist, abs_phi_longlist, abs_apar_longlist,
               abs_bpar_longlist, growth_rate_longlist, freq_longlist)

def examine_basic_gradient_scan():
    """ """
    sims_prefix = "sims/em_explicit_gradients_"

    outnc_longnames = [sims_prefix + "0" + ".out.nc",
                       sims_prefix + "1" + ".out.nc",
                       sims_prefix + "2" + ".out.nc",
                       sims_prefix + "3" + ".out.nc",
                       sims_prefix + "3_all_electron" + ".out.nc",
                       sims_prefix + "3_all_ion" + ".out.nc",
                       sims_prefix + "3_double_beta" + ".out.nc",
                      ]
    (gradient_longlist, z_longlist, abs_phi_longlist, abs_apar_longlist,
           abs_bpar_longlist, growth_rate_longlist, freq_longlist) = postprocess_folder_stella(outnc_longnames, make_plot=True)

    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # ax1.

    return

if __name__ == "__main__":
    print("Hello world")
    examine_basic_gradient_scan()
