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

def sort_with_idxs(input_list, idxs):
    """ """

    new_list = []
    for idx in idxs:
        new_list.append(input_list[idx])

    return new_list

def postprocess_folder_stella(outnc_longnames, param_scanned="gradients",
                            save_as_pickle=False, make_plot=False, sort=True):
    """For a folder, get Omega and beta, and save to pickle"""

    gradient_longlist = []
    growth_rate_longlist = []
    freq_longlist = []
    z_longlist = []
    abs_phi_longlist = []
    abs_apar_longlist = []
    abs_bpar_longlist = []
    time_longlist = []
    gammaom_longlist = []
    freqom_longlist = []

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
        max_abs_phi = np.max(abs_phi)
        abs_phi_longlist.append(abs_phi/max_abs_phi)
        abs_apar_longlist.append(abs_apar/max_abs_phi)
        abs_bpar_longlist.append(abs_bpar/max_abs_phi)
        time_longlist.append(time)
        gammaom_longlist.append(gammaom)
        freqom_longlist.append(freqom)
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
    if sort:
        sort_idxs = np.argsort(gradient_vals)
        gradient_vals = gradient_vals[sort_idxs]
        growth_rate_vals = growth_rate_vals[sort_idxs]
        freq_vals = freq_vals[sort_idxs]
        z_longlist = sort_with_idxs(z_longlist, sort_idxs)
        abs_phi_longlist = sort_with_idxs(abs_phi_longlist, sort_idxs)
        abs_apar_longlist = sort_with_idxs(abs_apar_longlist, sort_idxs)
        abs_bpar_longlist = sort_with_idxs(abs_bpar_longlist, sort_idxs)
        gammaom_longlist = sort_with_idxs(gammaom_longlist, sort_idxs)
        freqom_longlist = sort_with_idxs(freqom_longlist, sort_idxs)
        time_longlist = sort_with_idxs(time_longlist, sort_idxs)
    ## Save to a pickle
    if save_as_pickle:
        pickle_name = folder_longname + "/beta_gamma_omega.pickle"
        myfile = open(pickle_name, "wb")
        pickle.dump((gradient_vals, z_longlist, abs_phi_longlist, abs_apar_longlist,
               abs_bpar_longlist, growth_rate_vals, freq_vals,
               time_longlist, gammaom_longlist, freqom_longlist), myfile)
        myfile.close()

        return
    else:
        return (gradient_vals, z_longlist, abs_phi_longlist, abs_apar_longlist,
               abs_bpar_longlist, growth_rate_vals, freq_vals,
               time_longlist, gammaom_longlist, freqom_longlist)

def examine_basic_gradient_scan():
    """ """
    sims_prefix = "sims/em_explicit_gradients_"

    outnc_longnames = [sims_prefix + "0" + ".out.nc",       # 0
                       sims_prefix + "1" + ".out.nc",       # 1
                       sims_prefix + "2" + ".out.nc",       # 2
                       sims_prefix + "3" + ".out.nc",       # 3
                       sims_prefix + "3_all_electron" + ".out.nc",  # 4
                       sims_prefix + "3_all_ion" + ".out.nc",       # 5
                       sims_prefix + "3_double_beta" + ".out.nc",   # 6
                       "sims/em_implicit_gradients_3.out.nc"        # 7
                      ]
    (gradient_longlist, z_longlist, abs_phi_longlist, abs_apar_longlist,
           abs_bpar_longlist, growth_rate_longlist, freq_longlist,
           time_longlist, gammaom_longlist, freqom_longlist) = postprocess_folder_stella(outnc_longnames, make_plot=True, sort=False)

    fig = plt.figure()
    ax1 = fig.add_subplot(221) # omega
    ax2 = fig.add_subplot(223) # gamma
    ax3 = fig.add_subplot(322) # phi
    ax4 = fig.add_subplot(324) # apar
    ax5 = fig.add_subplot(326) # bpar

    my_linewidth = 3
    ylabel_fontsize = 20
    xlabel_fontsize = 20
    legend_fontsize = 12

    labels_list = [r"$a/L_{ni}=a/L_{ne}=a/L_{Ti}=a/L_{Te}=3, \beta=3\%$ (explicit)",
                   r"$a/L_{ni}=a/L_{Ti}=0, a/L_{ne}=a/L_{Te}=6, \beta=3\%$ (explicit)",
                   r"$a/L_{ni}=a/L_{Ti}=6, a/L_{ne}=a/L_{Te}=0, \beta=3\%$ (explicit)",
                   r"$a/L_{ni}=a/L_{ne}=a/L_{Ti}=a/L_{Te}=1.5, \beta=6\%$ (explicit)",
                   r"$a/L_{ni}=a/L_{ne}=a/L_{Ti}=a/L_{Te}=3, \beta=3\%$ (str_mirror implicit)",
                    ]
    sim_idxs = [3, 4, 5, 6, 7]

    # theta goes -pi to pi
    # each field period is 2/5 * pi
    # so nfield_period = 10 means zeta spans 4pi i.e. goes -2pi to 2pi
    # so nfield_period = 4 means zeta spans 8/5*pi i.e. goes -4/5*pi to 4/5*pi
    # so we need to go from ()
    explicit_zfactor = 4/(2*np.pi)
    implicit_zfactor = 1.6/(2*np.pi)
    conversion_factor = [explicit_zfactor, explicit_zfactor, explicit_zfactor,
                         explicit_zfactor, implicit_zfactor]
    linestyles=["-", "--", "-.", ":", "-"]
    my_linewidths=[5, 4, 4, 4, 2]
    for counter, sim_idx in enumerate(sim_idxs):
        z = z_longlist[sim_idx]
        abs_phi = abs_phi_longlist[sim_idx]
        abs_apar = abs_apar_longlist[sim_idx]
        abs_bpar = abs_bpar_longlist[sim_idx]
        time = time_longlist[sim_idx] ; time = time - time[-1]
        gammaom = gammaom_longlist[sim_idx]
        freqom = freqom_longlist[sim_idx]
        sim_label=labels_list[counter]

        ax1.plot(time, freqom, lw=my_linewidths[counter], label=sim_label, ls=linestyles[counter])
        ax2.plot(time, gammaom, lw=my_linewidths[counter], label=sim_label, ls=linestyles[counter])
        ax3.plot(z*conversion_factor[counter], abs_phi, lw=my_linewidths[counter], label=sim_label, ls=linestyles[counter])
        ax4.plot(z*conversion_factor[counter], abs_apar, lw=my_linewidths[counter], label=sim_label, ls=linestyles[counter])
        ax5.plot(z*conversion_factor[counter], abs_bpar, lw=my_linewidths[counter], label=sim_label, ls=linestyles[counter])

    ax2.set_ylim((0.5,4))
    ax1.set_ylim((-5, 5))
    ax1.set_ylabel(r"$\omega$", fontsize=ylabel_fontsize)
    ax2.set_ylabel(r"$\gamma$", fontsize=ylabel_fontsize)
    ax2.set_xlabel(r"$t-t_{end}$", fontsize=xlabel_fontsize)
    ax3.set_ylabel(r"$\vert \phi \vert$", fontsize=ylabel_fontsize)
    ax4.set_ylabel(r"$\vert A_\parallel \vert$", fontsize=ylabel_fontsize)
    ax5.set_ylabel(r"$\vert B_\parallel \vert$", fontsize=ylabel_fontsize)
    ax5.set_xlabel(r"$\zeta$", fontsize=xlabel_fontsize)
    ax2.legend(loc="best", fontsize=legend_fontsize)
    plt.tight_layout()
    plt.show()

    return

if __name__ == "__main__":
    print("Hello world")
    examine_basic_gradient_scan()
