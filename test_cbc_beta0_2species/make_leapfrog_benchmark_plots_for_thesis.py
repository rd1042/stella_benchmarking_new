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
stella_leapfrog_drifts_longname = "sims/stella_leapfrog_archer2/input"
gs2_basecase_longname = "sims/gs2_ky05/_0.0000"

def plot_phi_z_for_sim(ax1, theta, abs_phi, sim_label, sim_type="stella", plot_format=".eps",
                linewidth=1.0, linestyleoveride=False, color=False):
    """ """
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

    mylinewidths_1 = [10., 7, 5]
    mylinewidths_2 = [15., 9, 7]
    # mycols = ["black", "cyan", "red"]
    myaxlabelsize = 50.
    mylegendsize = 25
    mylabelsize = 40.
    mymarkersize = 20.
    my_ylabelsize = 30
    ## Plot of omega(t)
    fig1 = plt.figure(figsize=[12, 12])

    fig1_left = 0.16
    fig1_right = 0.98
    fig1_bottom = 0.14
    fig1_top = 0.98
    vspace = 0.05
    fig1_width = fig1_right - fig1_left
    fig1_height = (fig1_top - fig1_bottom - 2*vspace)/3
    row2_bottom = fig1_bottom + vspace + fig1_height
    row1_bottom = row2_bottom + vspace + fig1_height
    ax11 = fig1.add_axes((fig1_left, row1_bottom, fig1_width, fig1_height))
    ax12 = fig1.add_axes((fig1_left, row2_bottom, fig1_width, fig1_height))
    ax13 = fig1.add_axes((fig1_left, fig1_bottom, fig1_width, fig1_height))

    ## Plot of |phi|(t)
    fig2 = plt.figure(figsize=[12, 12])
    fig2_left = 0.17
    fig2_right = 0.985
    fig2_bottom = 0.12
    fig2_top = 0.97
    fig2_width = fig2_right - fig2_left
    fig2_height = fig2_top - fig2_bottom
    ax21 = fig2.add_axes((fig2_left, fig2_bottom, fig2_width, fig2_height))
    gamma_vals = []
    freq_vals = []

    gamma_llims = []
    gamma_ulims = []
    freq_llims = []
    freq_ulims = []

    linestyleoverides = ["-", "--", "-."]
    sim_data_pickle = "sim_summary.pickle"
    try:
        myfile = open(sim_data_pickle, "rb")
        [sim_data] = pickle.load(myfile)
        myfile.close()
        pickle_exists = True
    except FileNotFoundError:
        pickle_exists = False
        sim_data = []
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
        if pickle_exists:
            [time, freqom_final, gammaom_final, freqom, gammaom, gamma_stable, z, phi] = sim_data[sim_idx]
        else:
            time, freqom_final, gammaom_final, freqom, gammaom, gamma_stable = get_omega_data(sim_longname, "stella")
            z, real_phi, imag_phi = get_phiz_data(sim_longname, sim_type)
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
            phi = np.sqrt(real_phi*real_phi + imag_phi*imag_phi)

            ## Normalise s.t. max(abs_phi) = 1
            phi = phi/np.max(phi)
            sim_data.append([time, freqom_final, gammaom_final, freqom, gammaom, gamma_stable, z, phi])

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
        print("Omega = ", np.float(freqom_final), np.float(gamma_stable[-1]))
        time = np.array(time)
        freqom = np.array(freqom)
        gammaom = np.array(gammaom)
        gamma_stable = np.array(gamma_stable)
        ax11.plot(time, freqom, linewidth=mylinewidths_1[sim_idx], ls=linestyleoverides[sim_idx])
        ax12.plot(time, gammaom, linewidth=mylinewidths_1[sim_idx], ls=linestyleoverides[sim_idx])
        ax13.plot(time, gamma_stable, linewidth=mylinewidths_1[sim_idx], ls=linestyleoverides[sim_idx])
        plot_phi_z_for_sim(ax21, z, phi, sim_label, sim_type=sim_type, linewidth=mylinewidths_2[sim_idx],
                    linestyleoveride=linestyleoverides[sim_idx])#, color=mycols[sim_idx])


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
    # ax11.set_ylim(freq_llim, freq_ulim)
    # ax12.set_ylim(gamma_llim, gamma_ulim)
    ax11.set_ylim(0.16, 0.29)
    ax12.set_ylim(0.1, 0.25)
    ax13.set_ylim(0.165, 0.195)
    ax21.set_ylim(-0.05, 1.05)
    ax21.set_xlim(-6, 6)
    ax13.set_xlabel(r"$\tilde{t}$", fontsize=myaxlabelsize)
    ax11.set_ylabel(r"$\tilde{\omega}_{stella}$", fontsize=myaxlabelsize)
    ax12.set_ylabel(r"$\tilde{\gamma}_{stella}$", fontsize=myaxlabelsize)
    ax13.set_ylabel(r"$\tilde{\gamma}_{2}$", fontsize=myaxlabelsize)
    ax21.set_xlabel(r"$z/\pi$", fontsize=myaxlabelsize)
    ax21.set_ylabel(r"$\vert \tilde{\phi}_k \vert$", fontsize=myaxlabelsize)
    axes = [ax21]

    ax21.legend(loc="upper left", fontsize=mylegendsize, fancybox=True, framealpha=0.5)


    for ax in axes:
        ax.grid(True)
        #ax.legend(loc="best")
        ax.tick_params(axis='both', which='major', labelsize=mylabelsize, direction="in")
        ax.tick_params(axis='both', which='minor', direction="in")

    axes = [ax11, ax12, ax13]
    for ax in axes:
        ax.set_xlim(0,125)
        ax.grid(True)
        ax.tick_params(axis='both', which='major', labelsize=mylabelsize, direction="in")
        ax.tick_params(axis='both', which='minor', direction="in")
        ax.set_xticks([0, 40, 80, 120])

    ax13.set_xticklabels([r"$0$", r"$40$", r"$80$", r"$120$"], fontsize=mylabelsize )
    ax21.set_xticks([-6, -4, -2, 0, 2, 4, 6] )
    ax21.set_xticklabels([r"$-6$", r"$-4$", r"$-2$", r"$0$", r"$2$", r"$4$", r"$6$"], fontsize=mylabelsize )

    ax11.set_yticks([0.2, 0.25])
    ax11.set_yticklabels([r"$0.2$", r"$0.25$"], fontsize=my_ylabelsize )
    ax12.set_yticks([0.15, 0.2])
    ax12.set_yticklabels([r"$0.15$", r"$0.2$"], fontsize=my_ylabelsize )
    ax13.set_yticks([0.17, 0.18, 0.19])
    ax13.set_yticklabels([r"$0.17$", r"$0.18$", r"$0.19$"], fontsize=my_ylabelsize )
    for ax in [ax11, ax12]:
        ax.set_xticklabels([])

    #plt.show()
    fig1.savefig(save_name + "_omega" + ".eps")
    fig2.savefig(save_name + "_phi" + ".eps")
    plt.close(fig1)
    plt.close(fig2)

    if not pickle_exists:
        myfile = open(sim_data_pickle, "wb")
        pickle.dump([sim_data],myfile)
        myfile.close()
    return


def compare_stella_leapfrog_gs2_for_thesis():
    """ """

    make_comparison_plots_leapfrog_for_thesis([
            gs2_basecase_longname,
            stella_master_longname,
            stella_leapfrog_drifts_longname,
                    ],
                    [
                    r"GS2",#;$\Omega$=$0.202$+$0.187$i",
                    r"stella",#;$\Omega$=$0.206$+$0.184$i",
                    r"stella (Multipstep)"#";$\Omega$=$0.201$+$0.138$i",
                    # r"GS2;$\Omega$=$0.202$+$0.187$i",
                    # r"stella;$\Omega$=$0.206$+$0.184$i",
                    # r"stella (Multipstep);$\Omega$=$0.201$+$0.138$i",
                    ],
                    "./linear_leapfrog_benchmark",
                    sim_types=[
                    "gs2",
                    "stella",
                    "stella",
                    ],
                    )
    return

if __name__ == "__main__":
    compare_stella_leapfrog_gs2_for_thesis()
