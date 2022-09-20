"""Plot the linear test"""

import sys
sys.path.append("../postprocessing_tools")
from extract_sim_data import get_omega_data, get_phiz_data
import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
import sys
import pickle

stella_master_longname = "sims/stella_master_archer2/input"
stella_nperiod5_leapfrog_drifts_longname = "sims/stella_cmiller_es_2species_leapfrog_drifts/input5"
gs2_basecase_longname = "sims/gs2_ky05/_0.0000"

def plot_phi_z_for_sim(ax1, sim_longname, sim_label, sim_type="stella", plot_format=".eps",
                linewidth=1.0, linestyleoveride=False, color=False):
    """ """

    theta, real_phi, imag_phi = get_phiz_data(sim_longname, sim_type)
    ## Check values are finite
    if not(np.all(np.isfinite(real_phi)) and np.all(np.isfinite(imag_phi))):
        print("Error! phi contains non-finite values")
        print("sim_longname = ", sim_longname)
        return

    ## Combine real and imaginary parts to get abs_phi
    # If real and imaginary parts are large, it's possible that they'll
    # become non-finite when we square them. To avoid, perform some normalisation first
    normalisation = np.max(abs(real_phi))
    real_phi = real_phi/normalisation
    imag_phi = imag_phi/normalisation
    abs_phi = np.sqrt(real_phi*real_phi + imag_phi*imag_phi)

    ## Normalise s.t. max(abs_phi) = 1
    abs_phi = abs_phi/np.max(abs_phi)
    # Plot
    if not linestyleoveride:
        linestyle=next(linestyles2)
    else:
        linestyle = linestyleoveride
    if not color:
        ax1.plot(theta/np.pi, abs_phi, label=sim_label, ls=linestyle, lw=linewidth)
    else:
        ax1.plot(theta/np.pi, abs_phi, label=sim_label, ls=linestyle, lw=linewidth,c=color)

    return

def make_comparison_plots_leapfrog_for_thesis(sim_longnames, sim_labels, save_name, sim_types=[],
                          ):
    """Compare multiple simulations which have a single common input. Create the following
    plots:
    1) omega(t)
    2) Normalised |phi|(z)
    """

    mylinewidths = [15., 9, 7]
    mycols = ["black", "cyan", "red"]
    myaxlabelsize = 50.
    mylegendsize = 35
    mylabelsize = 40.
    mymarkersize = 20.

    ## Plot of omega(t)
    fig1 = plt.figure(figsize=[10, 12])
    ax11 = fig1.add_subplot(211)
    ax12 = fig1.add_subplot(212, sharex=ax11)

    ## Plot of |phi|(t)
    fig2 = plt.figure(figsize=[17, 12])
    ax21 = fig2.add_subplot(111)
    gamma_vals = []
    freq_vals = []

    gamma_llims = []
    gamma_ulims = []
    freq_llims = []
    freq_ulims = []

    linestyleoverides = ["-", "--", "-."]
    for sim_idx, sim_longname in enumerate(sim_longnames):
        sim_label = sim_labels[sim_idx]
        # Find out the sim types - if sim_types kwarg not specified,
        # assume all stella
        if len(sim_types) == 0:
            sim_type="stella"
        elif len(sim_types) == len(sim_longnames):
            sim_type = sim_types[sim_idx]
        else:
            print("Error! len(sim_longnames), len(sim_types) = ", len(sim_longnames), len(sim_types) )
            sys.exit()
        time, freqom_final, gammaom_final, freqom, gammaom, gamma_stable = get_omega_data(sim_longname, "stella")
        # print("time = ", np.array(time))
        # print("freqom_final = ", np.array(freqom_final))
        # print("gammaom_final = ", np.array(gammaom_final))
        # print("freqom = ", np.array(freqom))
        # print("gammaom = ", np.array(gammaom))
        # print("gamma_stable) = ", np.array(gamma_stable))
        gamma_llim = np.min(gammaom)
        gamma_ulim = np.max(gammaom)
        freq_llim = np.min(freqom)
        freq_ulim = np.max(freqom)
        print("Omega = ", np.float(freqom_final), np.float(gammaom_final))
        time = np.array(time)
        freqom = np.array(freqom)
        gammaom = np.array(gammaom)
        gamma_stable = np.array(gamma_stable)
        ax11.plot(time, freqom)
        ax12.plot(time, gammaom)
        plot_phi_z_for_sim(ax21, sim_longname, sim_label, sim_type=sim_type, linewidth=mylinewidths[sim_idx],
                    linestyleoveride=linestyleoverides[sim_idx], color=mycols[sim_idx])


        if (np.isfinite(gammaom_final) and np.isfinite(freqom_final)
                and np.isfinite(gamma_llim) and np.isfinite(gamma_ulim)
                and np.isfinite(freq_llim) and np.isfinite(freq_ulim) ):
            gamma_llims.append(gamma_llim)
            gamma_ulims.append(gamma_ulim)
            freq_llims.append(freq_llim)
            freq_ulims.append(freq_ulim)
            gamma_vals.append(gammaom_final)
            freq_vals.append(freqom_final)

    ## Set lims based on sim data
    gamma_llim = np.min(np.array(gamma_llims))/1.1
    gamma_ulim = np.max(np.array(gamma_ulims))*1.1
    freq_llim = np.min(np.array(freq_llims))/1.1
    freq_ulim = np.max(np.array(freq_ulims))*1.1
    ax11.set_ylim(freq_llim, freq_ulim)
    ax12.set_ylim(gamma_llim, gamma_ulim)
    ax21.set_ylim(-0.05, 1.19)
    ax21.set_xlim(-6, 6)
    ax12.set_xlabel(r"$t$", fontsize=myaxlabelsize)
    ax11.set_ylabel(r"$\omega$", fontsize=myaxlabelsize)
    ax12.set_ylabel(r"$\gamma$", fontsize=myaxlabelsize)
    ax21.set_xlabel(r"$z/\pi$", fontsize=myaxlabelsize)
    ax21.set_ylabel(r"$\vert \phi \vert$", fontsize=myaxlabelsize)
    axes = [ax21]

    ax21.legend(loc="upper left", fontsize=mylegendsize, fancybox=True, framealpha=0.5)


    for ax in axes:
        ax.grid(True)
        #ax.legend(loc="best")
        ax.tick_params(axis='both', which='major', labelsize=mylabelsize, direction="in")
        ax.tick_params(axis='both', which='minor', direction="in")

    for fig in [fig2]:
        fig.tight_layout()
    #plt.show()
    fig1.savefig(save_name + "_omega" + ".eps")
    fig2.savefig(save_name + "_phi" + ".eps")
    plt.close(fig1)
    plt.close(fig2)
    return


def compare_stella_leapfrog_gs2_for_thesis():
    """ """

    make_comparison_plots_leapfrog_for_thesis([
            stella_master_longname,
            # stella_nperiod5_leapfrog_drifts_longname,
            gs2_basecase_longname,
                    ],
                    [
                    r"stella;$\Omega$=$0.206$+$0.184$i",
                    # r"stella (Multipstep);$\Omega$=$0.12$+$0.14$i",
                    r"GS2;$\Omega$=$0.202$+$0.187$i",
                    ],
                    "./poster_leapfrog",
                    sim_types=[
                    "stella",
                    # "stella",
                    "gs2",
                    ],
                    )
    return

if __name__ == "__main__":
    compare_stella_leapfrog_gs2_for_thesis()
