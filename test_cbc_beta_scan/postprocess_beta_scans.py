"""Postprocessing simulations for GS2 and stella, so we can compare linear simulations
(e.g. beta scans).
For each simulation, we want to extrac omega(t=tfinal)
It would also be handy to plot Omega(t), to check convergence
And fields(t)
 """

import sys
import pickle
import numpy as np
import glob
import re
import matplotlib.pyplot as plt

sys.path.append("../postprocessing_tools")
import helper_linear_sims as help_lin
from helper_ncdf import view_ncdf_variables
from extract_sim_data import get_omega_data, get_phiz_data, get_aparz_data, get_bparz_data

def get_beta_from_outnc_longname(outnc_longname):
    """ """

    sim_shortname = re.split("/", outnc_longname)[-1]
    sim_shortname = re.split(".out.nc", sim_shortname)[0]
    beta_str = re.split("beta_", sim_shortname)[-1]
    return float(beta_str)

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

def postprocess_folder_gs2(folder_shortname, param_scanned="beta"):
    """For a folder, get Omega and beta, and save to pickle"""

    folder_longname = "sims/" + folder_shortname
    if param_scanned != "beta":
        print("Error! Only works for param_scanned=beta")
        sys.exit()

    beta_longlist = []
    growth_rate_longlist = []
    freq_longlist = []

    outnc_longnames = glob.glob(folder_longname + "/*.out.nc")

    for outnc_longname in outnc_longnames:

        sim_shortname = re.split("/", outnc_longname)[-1]
        sim_shortname = re.split(".out.nc", sim_shortname)[0]
        ## Get beta
        beta_longlist.append(get_beta_from_outnc_longname(outnc_longname))
        (fitted_growth_rate, growth_rate_error,
         converged_frequency, freq_error) = help_lin.calculate_omega_and_plot_for_single(outnc_longname,
                        save_path=(folder_longname + "/"),
                        plot_growth_rates=True,
                        figtitle=sim_shortname)

        growth_rate_longlist.append(fitted_growth_rate)
        freq_longlist.append(converged_frequency)

    ## Sort beta and Omega
    beta_vals = np.array(beta_longlist)
    growth_rate_vals = np.array(growth_rate_longlist)
    freq_vals = np.array(freq_longlist)
    sort_idxs = np.argsort(beta_vals)
    beta_vals = beta_vals[sort_idxs]
    growth_rate_vals = growth_rate_vals[sort_idxs]
    freq_vals = freq_vals[sort_idxs]

    ## Save to a pickle
    pickle_name = folder_longname + "/beta_gamma_omega.pickle"
    myfile = open(pickle_name, "wb")
    pickle.dump([beta_vals, growth_rate_vals, freq_vals], myfile)
    myfile.close()
    return

def postprocess_folder_stella(folder_shortname, param_scanned="beta"):
    """For a folder, get Omega and beta, and save to pickle"""

    folder_longname = "sims/" + folder_shortname
    if param_scanned != "beta":
        print("Error! Only works for param_scanned=beta")
        sys.exit()

    beta_longlist = []
    growth_rate_longlist = []
    freq_longlist = []

    ## stella simulations tend to be expensive, so might timeout. Can we
    ## still get mode structures, omega(t) from the .out.nc?
    ## partially-completed .out.nc files contain:
    ## ['code_info', 'nproc', 'nmesh', 'ntubes', 'nkx', 'nky', 'nzed_tot',
    ## 'nspecies', 'nmu', 'nvpa_tot', 't', 'charge', 'mass', 'dens', 'temp',
    ## 'tprim', 'fprim', 'vnew', 'type_of_species', 'theta0', 'kx', 'ky', 'mu',
    ## 'vpa', 'zed', 'bmag', 'gradpar', 'gbdrift', 'gbdrift0', 'cvdrift',
    ## 'cvdrift0', 'kperp2', 'gds2', 'gds21', 'gds22', 'grho', 'jacob', 'q',
    ## 'beta', 'shat', 'jtwist', 'drhodpsi', 'phi2', 'bpar2', 'input_file']
    outnc_longnames = glob.glob(folder_longname + "/*.out.nc")

    for outnc_longname in outnc_longnames:

        sim_shortname = re.split("/", outnc_longname)[-1]
        sim_shortname = re.split(".out.nc", sim_shortname)[0]
        sim_longname = re.split(".out.nc", outnc_longname)[0]
        # view_ncdf_variables(outnc_longname)
        ## Get beta
        beta_longlist.append(get_beta_from_outnc_longname(outnc_longname))
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
            abs_bpar = help_lin.get_abs(real_bpar, imag_apar)
        except None:
            print("None")
        make_omega_time_plot_for_stella(sim_longname, time, freqom, gammaom, gamma_stable)
        make_field_z_plot_for_stella(sim_longname, z, abs_phi, abs_apar, abs_bpar)
        # (fitted_growth_rate, growth_rate_error,
        #  converged_frequency, freq_error) = help_lin.calculate_omega_and_plot_for_single(outnc_longname,
        #                 save_path=(folder_longname + "/"),
        #                 plot_growth_rates=True,
        #                 figtitle=sim_shortname)

        growth_rate_longlist.append(gamma_stable[-1])
        freq_longlist.append(freqom_final)

    ## Sort beta and Omega
    beta_vals = np.array(beta_longlist)
    growth_rate_vals = np.array(growth_rate_longlist)
    freq_vals = np.array(freq_longlist)
    sort_idxs = np.argsort(beta_vals)
    beta_vals = beta_vals[sort_idxs]
    growth_rate_vals = growth_rate_vals[sort_idxs]
    freq_vals = freq_vals[sort_idxs]

    ## Save to a pickle
    pickle_name = folder_longname + "/beta_gamma_omega.pickle"
    myfile = open(pickle_name, "wb")
    pickle.dump([beta_vals, growth_rate_vals, freq_vals], myfile)
    myfile.close()
    return

def postprocess_folder(folder_shortname, sim_type, param_scanned="beta"):
    """ """
    if sim_type == "gs2":
        postprocess_folder_gs2(folder_shortname, param_scanned=param_scanned)
    elif sim_type == "stella":
        postprocess_folder_stella(folder_shortname, param_scanned=param_scanned)

    return

def make_plot_for_single_sim(sim_longname):
    """ """
    sim_shortname = re.split("/", sim_longname)[-1]
    # view_ncdf_variables(outnc_longname)
    ## Get beta
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
        abs_bpar = help_lin.get_abs(real_bpar, imag_apar)
    except None:
        print("None")
    make_omega_time_plot_for_stella(sim_longname, time, freqom, gammaom, gamma_stable)
    make_field_z_plot_for_stella(sim_longname, z, abs_phi, abs_apar, abs_bpar)
    return

if __name__ == "__main__":
    # postprocess_folder("gs2_beta_scan_ky_05_np4_nt64_ng12_ne24_fapar1_fbpar1", "gs2")
    # postprocess_folder("gs2_beta_scan_ky_05_np3_nt48_ng8_ne18_fapar1_fbpar1", "gs2")
    postprocess_folder("stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1", "stella")
    # postprocess_folder("stella_beta_scan_ky_05_np2_nt64_nvpa24_nmu18_fapar1_fbpar1", "stella")
    # postprocess_folder("stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1", "stella")
    # postprocess_folder("stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1", "stella")
    # compare_em_old = "old_em_stella_comparison/electromagnetic/"
    # postprocess_folder(compare_em_old + "stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2", "stella")
    # postprocess_folder(compare_em_old + "stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2", "stella")
    # postprocess_folder(compare_em_old + "stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2", "stella")
    # compare_em_new = "old_em_stella_comparison/electromagnetic-new/"
    # postprocess_folder(compare_em_new + "stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2", "stella")
    # postprocess_folder(compare_em_new + "stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2", "stella")
    # postprocess_folder(compare_em_new + "stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2", "stella")
    # compare_em_new_impl_str = "old_em_stella_comparison/electromagnetic-new_stream_implicit/"
    # postprocess_folder(compare_em_new_impl_str + "stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2", "stella")
    # postprocess_folder(compare_em_new_impl_str + "stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2", "stella")
    # postprocess_folder(compare_em_new_impl_str + "stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2", "stella")
    #
    # postprocess_folder("stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_streaming_implicit", "stella")
    # postprocess_folder("stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_streaming_implicit", "stella")
    # postprocess_folder("stella_beta_scan_ky_05_np2_nt64_nvpa24_nmu18_fapar1_fbpar1_streaming_implicit", "stella")
    # postprocess_folder("stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1_streaming_implicit", "stella")
    # gfort_compare_em_new = "viking_gfortran_build/old_em_stella_comparison/electromagnetic-new/"
    # postprocess_folder(gfort_compare_em_new + "stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2", "stella")
    # gfort_compare_em_new_impl_str = "viking_gfortran_build/old_em_stella_comparison/electromagnetic-new_stream_implicit/"
    # postprocess_folder(gfort_compare_em_new_impl_str + "stella_beta_scan_ky_05_np2_nt128_nvpa18_nmu12_fapar1_fbpar1_2", "stella")
    # postprocess_folder(gfort_compare_em_new_impl_str + "stella_beta_scan_ky_05_np2_nt256_nvpa18_nmu12_fapar1_fbpar1_2", "stella")
    # postprocess_folder(gfort_compare_em_new_impl_str + "stella_beta_scan_ky_05_np2_nt512_nvpa18_nmu12_fapar1_fbpar1_2", "stella")
    # gfort_build = "viking_gfortran_build/"
    # postprocess_folder(gfort_build + "stella_beta_scan_ky_05_np2_nt32_nvpa18_nmu12_fapar1_fbpar1_streaming_implicit", "stella")
    # postprocess_folder(gfort_build + "stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_streaming_implicit", "stella")
    # make_plot_for_single_sim("sims/viking_gfortran_build/stella_beta_scan_ky_05_np2_nt64_nvpa18_nmu12_fapar1_fbpar1_streaming_implicit/" +
    #                          "beta_0.04000_investigation")
    # postprocess_folder(gfort_build + "stella_beta_scan_ky_05_np2_nt64_nvpa24_nmu18_fapar1_fbpar1_streaming_implicit", "stella")
    # postprocess_folder(gfort_build + "stella_beta_scan_ky_05_np4_nt64_nvpa36_nmu24_fapar1_fbpar1_streaming_implicit", "stella")
