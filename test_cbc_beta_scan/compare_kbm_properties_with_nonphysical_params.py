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

sys.path.append("../postprocessing_tools")
import helper_linear_sims as help_lin
from helper_ncdf import view_ncdf_variables
from extract_sim_data import get_omega_data, get_phiz_data, get_aparz_data, get_bparz_data
from plotting_helper import plot_phi_z_for_sim, plot_apar_z_for_sim, plot_bpar_z_for_sim
comparison_folder = "sims/beta_0.04000_investigation/"

def compare_all_sims():
    """ """
    outnc_longnames = glob.glob(comparison_folder + "*.out.nc")
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
            else:
                sim_type = "stella"

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
    for sim_idx, sim_shortname in enumerate(sim_shortnames):
        ax1.plot(tvals_list[sim_idx]-tend_list[sim_idx], omega_list[sim_idx], label=sim_shortname)
        ax2.plot(tvals_list[sim_idx]-tend_list[sim_idx], gamma_list[sim_idx], label=sim_shortname)

    max_omega_final = np.max(np.array(omega_final_list))
    min_omega_final = np.min(np.array(omega_final_list))
    max_gamma_final = np.max(np.array(gamma_final_list))
    min_gamma_final = np.min(np.array(gamma_final_list))
    gamma_buffer = 0.01
    omega_buffer = 0.01
    ax1.legend(loc="best")
    ax2.set_xlabel(r"$t$")
    ax1.set_ylabel(r"$\omega$")
    ax2.set_ylabel(r"$\gamma$")
    ax1.set_ylim(min_omega_final-omega_buffer, max_omega_final+omega_buffer)
    ax2.set_ylim(min_gamma_final-gamma_buffer, max_gamma_final+gamma_buffer)

    for ax in [ax1, ax2]:
        ax.grid(True)

    plt.show()

if __name__ == "__main__":
    print("Hello world")
    compare_all_sims()
