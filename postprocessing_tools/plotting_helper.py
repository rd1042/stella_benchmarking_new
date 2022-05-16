""" """

from extract_sim_data import get_omega_data, get_phiz_data, get_aparz_data, get_bparz_data
from extract_sim_data import find_apar_phi_ratio, find_bpar_phi_ratio
from extract_sim_data import get_gs2_omega_from_plaintext, get_omega_data_gs2_outnc
from helper_ncdf import view_ncdf_variables, extract_data_from_ncdf, extract_data_from_ncdf_with_xarray
import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
import sys
import pickle

## Define linestyles
linestyles1=cycle(["-", "--", "-.", ":", (0, (3, 1, 1, 1, 1, 1))])
linestyles2=cycle(["-", "--", "-.", ":", (0, (3, 1, 1, 1, 1, 1))])


def plot_omega_t_for_sim(ax1, ax2, sim_longname, sim_label, sim_type="stella"):
    """ """
    print("sim_longname = ", sim_longname)
    time, freqom_final, gammaom_final, freqom, gammaom, gamma_stable = get_omega_data(sim_longname, sim_type)
    half_len = int(len(time)/2)
    freqom_half = freqom[half_len:]

    gammaom_half = gammaom[half_len:]
    # gammaom_half = gamma_stable[half_len:]
    # Get an estimate for the ylims by looking at the max/min of the second half
    # of the frequencuy and gamma data
    if np.isfinite(freqom_final) and np.isfinite(gammaom_final):
        linestyle=next(linestyles1)
        gamma_llim = (np.min(gammaom_half)); gamma_ulim = (np.max(gammaom_half))
        freq_llim = (np.min(freqom_half)); freq_ulim = np.max(freqom_half)
        ax1.plot(time, freqom, label=sim_label, ls=linestyle)
        ax2.plot(time, gammaom, label=sim_label, ls=linestyle)
        #ax2.plot(time, gamma_stable, label=sim_label + " gamma_stable")
    else:
        gamma_llim = np.NaN; gamma_ulim = np.NaN; freq_llim = np.NaN; freq_ulim = np.NaN

    return gammaom_final, freqom_final, gamma_llim, gamma_ulim, freq_llim, freq_ulim

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

def plot_apar_z_for_sim(ax1, sim_longname, sim_label, sim_type="stella", linewidth=1.0,
             linestyleoveride=False):
    """ """

    theta, real_apar, imag_apar = get_aparz_data(sim_longname, sim_type)
    ## Check values are finite
    if not(np.all(np.isfinite(real_apar)) and np.all(np.isfinite(imag_apar))):
        print("Error! apar contains non-finite values")
        print("sim_longname = ", sim_longname)
        return

    ## Combine real and imaginary parts to get abs_apar
    # If real and imaginary parts are large, it's possible that they'll
    # become non-finite when we square them. To avoid, perform some normalisation first
    normalisation = np.max(abs(real_apar))
    real_apar = real_apar/normalisation
    imag_apar = imag_apar/normalisation
    abs_apar = np.sqrt(real_apar*real_apar + imag_apar*imag_apar)

    ## Normalise s.t. max(abs_apar) = 1
    abs_apar = abs_apar/np.max(abs_apar)
    # Plot
    if not linestyleoveride:
        linestyle=next(linestyles2)
    else:
        linestyle = linestyleoveride
    ax1.plot(theta/np.pi, abs_apar, label=sim_label, ls=linestyle, lw=linewidth)

    return

def plot_bpar_z_for_sim(ax1, sim_longname, sim_label, sim_type="stella"):
    """ """

    theta, real_bpar, imag_bpar = get_bparz_data(sim_longname, sim_type)
    ## Check values are finite
    if not(np.all(np.isfinite(real_bpar)) and np.all(np.isfinite(imag_bpar))):
        print("Error! bpar contains non-finite values")
        print("sim_longname = ", sim_longname)
        return

    ## Combine real and imaginary parts to get abs_bpar
    # If real and imaginary parts are large, it's possible that they'll
    # become non-finite when we square them. To avoid, perform some normalisation first
    normalisation = np.max(abs(real_bpar))
    real_bpar = real_bpar/normalisation
    imag_bpar = imag_bpar/normalisation
    abs_bpar = np.sqrt(real_bpar*real_bpar + imag_bpar*imag_bpar)

    ## Normalise s.t. max(abs_bpar) = 1
    abs_bpar = abs_bpar/np.max(abs_bpar)
    # Plot
    linestyle=next(linestyles2)

    ax1.plot(theta/np.pi, abs_bpar, label=sim_label, ls=linestyle)

    return

def make_comparison_plots(sim_longnames, sim_labels, save_name, sim_types=[],
                          plot_apar=False, plot_bpar=False, plot_format=".eps", show_fig=False):
    """Compare multiple simulations which have a single common input. Create the following
    plots:
    1) omega(t)
    2) Normalised |phi|(z)
    """
    print("sim_longnames = ", sim_longnames)
    ## Plot of omega(t)
    fig1 = plt.figure(figsize=[10, 12])
    ax11 = fig1.add_subplot(211)
    ax12 = fig1.add_subplot(212, sharex=ax11)

    ## Plot of |phi|(t)
    fig2 = plt.figure(figsize=[8, 8])
    ax21 = fig2.add_subplot(111)
    if plot_apar:
        fig3 = plt.figure(figsize=[8, 8])
        ax31 = fig3.add_subplot(111)
    if plot_bpar:
        fig4 = plt.figure(figsize=[8, 8])
        ax41 = fig4.add_subplot(111)
    gamma_vals = []
    freq_vals = []

    gamma_llims = []
    gamma_ulims = []
    freq_llims = []
    freq_ulims = []

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
        gammaom_final, freqom_final, gamma_llim, gamma_ulim, \
            freq_llim, freq_ulim = plot_omega_t_for_sim(ax11, ax12, sim_longname, sim_label, sim_type=sim_type)
        plot_phi_z_for_sim(ax21, sim_longname, sim_label, sim_type=sim_type)
        if plot_apar:
            apar_ratio = find_apar_phi_ratio(sim_longname, sim_type)
            print("sim_longname, apar ratio = ", sim_longname, apar_ratio)
            plot_apar_z_for_sim(ax31, sim_longname, sim_label, sim_type=sim_type)
        if plot_bpar:
            plot_bpar_z_for_sim(ax41, sim_longname, sim_label, sim_type=sim_type)
            bpar_ratio = find_bpar_phi_ratio(sim_longname, sim_type)
            print("sim_longname, bpar ratio = ", sim_longname, bpar_ratio)

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
    ax21.set_ylim(-0.05, 1.05)
    ax12.set_xlabel(r"$t$")
    ax11.set_ylabel(r"$\omega$")
    ax12.set_ylabel(r"$\gamma$")
    ax21.set_xlabel(r"$\theta/\pi$")
    ax21.set_ylabel(r"$\vert \phi \vert$")
    axes = [ax11, ax12, ax21]
    if plot_apar:
        axes.append(ax31)
        ax31.set_ylim(-0.05, 1.05)
        ax31.set_xlabel(r"$\theta/\pi$")
        ax31.set_ylabel(r"$\vert A_\parallel \vert$")
    if plot_bpar:
        axes.append(ax41)
        ax41.set_ylim(-0.05, 1.05)
        ax41.set_xlabel(r"$\theta/\pi$")
        ax41.set_ylabel(r"$\vert B_\parallel \vert$")
    for ax in axes:
        ax.grid(True)
        ax.legend(loc="best")
    if show_fig:
        plt.show()
    else:
        fig1.savefig(save_name + "_omega" + plot_format)
        fig2.savefig(save_name + "_phi" + plot_format)
        plt.close(fig1)
        plt.close(fig2)
    if plot_apar:
        if not show_fig:
            fig3.savefig(save_name + "_apar" + plot_format)
            plt.close(fig3)
    if plot_bpar:
        if not show_fig:
            fig4.savefig(save_name + "_bpar" + plot_format)
            plt.close(fig4)
    return

def make_comparison_plots_for_poster(sim_longnames, sim_labels, save_name, sim_types=[],
                          ):
    """Compare multiple simulations which have a single common input. Create the following
    plots:
    1) omega(t)
    2) Normalised |phi|(z)
    """

    mylinewidth = 6.
    myaxlabelsize = 50.
    mylegendsize = 35.
    mylabelsize = 40.
    mymarkersize = 20.

    ## Plot of omega(t)
    fig1 = plt.figure(figsize=[10, 12])
    ax11 = fig1.add_subplot(211)
    ax12 = fig1.add_subplot(212, sharex=ax11)

    ## Plot of |phi|(t)
    fig2 = plt.figure(figsize=[8, 8])
    ax21 = fig2.add_subplot(111)
    fig3 = plt.figure(figsize=[8, 8])
    ax31 = fig3.add_subplot(111)
    gamma_vals = []
    freq_vals = []

    gamma_llims = []
    gamma_ulims = []
    freq_llims = []
    freq_ulims = []

    [stella_longname, gs2_longname] = sim_longnames
    linestyleoverides = ["-", "-."]
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
        gammaom_final, freqom_final, gamma_llim, gamma_ulim, \
            freq_llim, freq_ulim = plot_omega_t_for_sim(ax11, ax12, sim_longname, sim_label, sim_type=sim_type)
        plot_phi_z_for_sim(ax21, sim_longname, sim_label, sim_type=sim_type, linewidth=mylinewidth, linestyleoveride=linestyleoverides[sim_idx])
        apar_ratio = find_apar_phi_ratio(sim_longname, sim_type)
        print("sim_longname, apar ratio = ", sim_longname, apar_ratio)
        plot_apar_z_for_sim(ax31, sim_longname, sim_label, sim_type=sim_type, linewidth=mylinewidth, linestyleoveride=linestyleoverides[sim_idx])


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
    ax21.set_ylim(-0.05, 1.05)
    ax12.set_xlabel(r"$t$", fontsize=myaxlabelsize)
    ax11.set_ylabel(r"$\omega$", fontsize=myaxlabelsize)
    ax12.set_ylabel(r"$\gamma$", fontsize=myaxlabelsize)
    ax21.set_xlabel(r"$\theta/\pi$", fontsize=myaxlabelsize)
    ax21.set_ylabel(r"$\vert \phi \vert$", fontsize=myaxlabelsize)
    axes = [ax11, ax12, ax21]
    axes.append(ax31)
    ax31.set_ylim(-0.05, 1.05)
    ax31.set_xlabel(r"$\theta/\pi$", fontsize=myaxlabelsize)
    ax31.set_ylabel(r"$\vert A_\parallel \vert$", fontsize=myaxlabelsize)

    #ax21.legend(loc="upper right", fontsize=mylegendsize)
    #ax21.legend(loc="upper left", fontsize=mylegendsize, fancybox=True, framealpha=0.5)

    for ax in axes:
        ax.grid(True)
        ax.set_xlim(-6, 6)
        #ax.legend(loc="best")
        ax.tick_params(axis='both', which='major', labelsize=mylabelsize, direction="in")
        ax.tick_params(axis='both', which='minor', direction="in")

    for fig in [fig1, fig2, fig3]:
        fig.tight_layout()

    fig1.savefig(save_name + "_omega" + ".png")
    fig2.savefig(save_name + "_phi" + ".png")
    plt.close(fig1)
    plt.close(fig2)
    fig3.savefig(save_name + "_apar" + ".png")
    plt.close(fig3)
    return

def make_comparison_plots_leapfrog_poster(sim_longnames, sim_labels, save_name, sim_types=[],
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
    #fig1 = plt.figure(figsize=[10, 12])
    # ax11 = fig1.add_subplot(211)
    # ax12 = fig1.add_subplot(212, sharex=ax11)

    ## Plot of |phi|(t)
    fig2 = plt.figure(figsize=[17, 12])
    ax21 = fig2.add_subplot(111)
    gamma_vals = []
    freq_vals = []

    # gamma_llims = []
    # gamma_ulims = []
    # freq_llims = []
    # freq_ulims = []

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
        # gammaom_final, freqom_final, gamma_llim, gamma_ulim, \
        #     freq_llim, freq_ulim = plot_omega_t_for_sim(ax11, ax12, sim_longname, sim_label, sim_type=sim_type)
        #print("Omega = ", freqom_final, gammaom_final )
        plot_phi_z_for_sim(ax21, sim_longname, sim_label, sim_type=sim_type, linewidth=mylinewidths[sim_idx],
                    linestyleoveride=linestyleoverides[sim_idx], color=mycols[sim_idx])


        # if (np.isfinite(gammaom_final) and np.isfinite(freqom_final)
        #         and np.isfinite(gamma_llim) and np.isfinite(gamma_ulim)
        #         and np.isfinite(freq_llim) and np.isfinite(freq_ulim) ):
        #     gamma_llims.append(gamma_llim)
        #     gamma_ulims.append(gamma_ulim)
        #     freq_llims.append(freq_llim)
        #     freq_ulims.append(freq_ulim)
        #     gamma_vals.append(gammaom_final)
        #     freq_vals.append(freqom_final)

    ## Set lims based on sim data
    # gamma_llim = np.min(np.array(gamma_llims))/1.1
    # gamma_ulim = np.max(np.array(gamma_ulims))*1.1
    # freq_llim = np.min(np.array(freq_llims))/1.1
    # freq_ulim = np.max(np.array(freq_ulims))*1.1
    # ax11.set_ylim(freq_llim, freq_ulim)
    # ax12.set_ylim(gamma_llim, gamma_ulim)
    ax21.set_ylim(-0.05, 1.19)
    ax21.set_xlim(-6, 6)
    # ax12.set_xlabel(r"$t$", fontsize=myaxlabelsize)
    # ax11.set_ylabel(r"$\omega$", fontsize=myaxlabelsize)
    # ax12.set_ylabel(r"$\gamma$", fontsize=myaxlabelsize)
    ax21.set_xlabel(r"$\theta/\pi$", fontsize=myaxlabelsize)
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
    #fig1.savefig(save_name + "_omega" + ".eps")
    fig2.savefig(save_name + "_phi" + ".png")
    #plt.close(fig1)
    plt.close(fig2)
    return


def calculate_omega_for_param_scan(folder_name, param, param_key, file_regex=False,
                                   path=None, use_lin_filepath=False, **kwargs):
    """Extract output files from selected folder. For each output file, find the
    value of the scanned parameter and calculate Omega."""
    print("folder_name = ", folder_name)
    # print("Check your slashes")
    # sys.exit()
    # Extract the files
    if path is None:
        if use_lin_filepath:
            path = LIN_SIM_PATH + folder_name + "/"
        elif file_regex == False:
            path = folder_name + '/'
        else:
            path = folder_name + '/' + file_regex
    output_file_list = glob.glob(path + "*.out.nc")
    param_list = []
    frequency_list = []
    growth_rate_list = []
    freq_errors = []
    growth_rate_errors = []
    if len(output_file_list) == 0:
        print("output_file_list = ", output_file_list)
        sys.exit("Empty output file list, aborting")
    for i, output_file in enumerate(output_file_list):

        # Extract the name of the simulation from its full path
        m = re.search(path +'(.+?)\.out\.nc', output_file)
        found = m.group(1)
        save_loc = re.split((found+".out.nc"), output_file)[0]
        print("save_loc = ", save_loc)
        fitted_growth_rate, growth_rate_error, converged_frequency, freq_error = calculate_omega_and_plot_for_single(output_file,
                        figtitle=found, folder_name=folder_name, save_path=save_loc, **kwargs)
        frequency_list.append(converged_frequency)
        freq_errors.append(freq_error)
        growth_rate_list.append(fitted_growth_rate)
        growth_rate_errors.append(growth_rate_error)
        print("output_file, fitted_growth_rate = ", output_file, fitted_growth_rate)
        try:    # Bob: should fix this
            param_list.append(help_ncdf.extract_data_from_ncdf(output_file, param_key))
        except KeyError:
            print("Error!")
            sys.exit()
    return (param_list, frequency_list, growth_rate_list,
            freq_errors, growth_rate_errors)


def make_beta_scan_plots(stella_longnames, gs2_longnames, beta_vals, save_name,
                         gs2_pickle=None,  plot_apar=False, plot_bpar=False, plot_format=".eps"):
    """Construct beta scans in the case where we have a set of stella sims and a
    set of GS2 sims."""

    ## Plot of omega(beta)
    fig2 = plt.figure(figsize=[10, 12])
    ax21 = fig2.add_subplot(211)
    ax22 = fig2.add_subplot(212, sharex=ax21)
    gamma_vals = []
    freq_vals = []
    final_beta_vals = []

    stella_gamma_vals = []
    stella_freq_vals = []
    gs2_gamma_vals = []
    gs2_freq_vals = []

    for sim_idx, stella_longname in enumerate(stella_longnames):
        ## Plot of omega(t)
        if len(gs2_longnames) == len(stella_longnames):
            make_comparison_plots([stella_longname, gs2_longnames[sim_idx]],
                                  ["stella beta = " + str(beta_vals[sim_idx]),
                                  "GS2 beta = " + str(beta_vals[sim_idx])
                                  ],
                                  (save_name + "_beta = " + str(beta_vals[sim_idx])),
                                  sim_types=["stella", "gs2"],
                                  plot_apar=plot_apar, plot_bpar=plot_bpar, plot_format=plot_format)

        time, freqom_final, gammaom_final, freqom, gammaom, gamma_stable = get_omega_data(stella_longname, "stella")
        if (np.isfinite(gammaom_final) and np.isfinite(freqom_final)):
            stella_gamma_vals.append(gammaom_final)
            stella_freq_vals.append(freqom_final)
            final_beta_vals.append(beta_vals[sim_idx])
        if len(gs2_longnames) == len(stella_longnames):
            try:
                freqom_final_gs2, gammaom_final_gs2 = get_gs2_omega_from_plaintext(gs2_longnames[sim_idx])
            except FileNotFoundError:
                t, freqom_final_gs2, gammaom_final_gs2, freq, gam, gam_stable = get_omega_data_gs2_outnc(gs2_longnames[sim_idx])


            if (np.isfinite(freqom_final_gs2) and np.isfinite(gammaom_final_gs2)):
                gs2_gamma_vals.append(gammaom_final_gs2)
                gs2_freq_vals.append(freqom_final_gs2)

    ax21.plot(final_beta_vals, stella_freq_vals, label="stella")
    ax22.plot(final_beta_vals, stella_gamma_vals, label="stella")
    ax21.scatter(final_beta_vals, stella_freq_vals, c="black", s=15., marker="x")
    ax22.scatter(final_beta_vals, stella_gamma_vals, c="black", s=15., marker="x")
    if len(gs2_longnames) == len(stella_longnames):
        ax21.plot(final_beta_vals, gs2_freq_vals, label="GS2")
        ax22.plot(final_beta_vals, gs2_gamma_vals, label="GS2")
        ax21.scatter(final_beta_vals, gs2_freq_vals, c="red", s=15., marker="x")
        ax22.scatter(final_beta_vals, gs2_gamma_vals, c="red", s=15., marker="x")

    if gs2_pickle is not None:
        myfile = open(gs2_pickle, 'rb')
        gs2_dict = pickle.load(myfile)  # contains 'beta', 'frequency', 'growth rate'
        myfile.close()
        gs2_beta = gs2_dict['beta']; gs2_frequency=gs2_dict['frequency']; gs2_growth_rate = gs2_dict['growth rate']
        # Sort the beta vals into ascending order
        gs2_beta = np.array(gs2_beta).flatten()
        sort_idxs = np.argsort(gs2_beta)
        gs2_beta = np.array(gs2_beta)[sort_idxs]
        gs2_frequency = np.array(gs2_frequency)[sort_idxs]
        gs2_growth_rate = np.array(gs2_growth_rate)[sort_idxs]
        ax21.plot(gs2_beta, gs2_frequency, label="GS2")
        ax22.plot(gs2_beta, gs2_growth_rate, label="GS2")

    ax21.set_ylabel(r"$\omega$")
    ax22.set_ylabel(r"$\gamma$")
    ax22.set_xlabel(r"$\beta$")

    for ax in [ax21, ax22]:
        ax.grid(True)
        ax.legend(loc="best")

    fig2.savefig(save_name + "_omega_beta" + plot_format)
    plt.close(fig2)

    return


def make_beta_scan_plots_for_poster(stella_longnames, beta_vals, save_name,
                         gs2_pickle=None):
    """Construct beta scans in the case where we have a set of stella sims and a
    set of GS2 sims."""
    mylinewidth = 6.
    myaxlabelsize = 60.
    mylegendsize = 40.
    mylabelsize = 40.
    mymarkersize = 20.

    ## Plot of omega(beta)
    fig2 = plt.figure(figsize=[14, 12])
    ax21 = fig2.add_subplot(211)
    ax22 = fig2.add_subplot(212, sharex=ax21)
    gamma_vals = []
    freq_vals = []
    final_beta_vals = []

    stella_gamma_vals = []
    stella_freq_vals = []
    gs2_gamma_vals = []
    gs2_freq_vals = []

    for sim_idx, stella_longname in enumerate(stella_longnames):
        ## Plot of omega(t)

        time, freqom_final, gammaom_final, freqom, gammaom, gamma_stable = get_omega_data(stella_longname, "stella")
        if (np.isfinite(gammaom_final) and np.isfinite(freqom_final)):
            stella_gamma_vals.append(gammaom_final)
            stella_freq_vals.append(freqom_final)
            final_beta_vals.append(beta_vals[sim_idx])

    ax21.plot(final_beta_vals, stella_freq_vals, label="stella", lw=mylinewidth)
    ax22.plot(final_beta_vals, stella_gamma_vals, label="stella", lw=mylinewidth)
    # ax21.scatter(final_beta_vals, stella_freq_vals, c="black", s=mymarkersize, marker="x")
    # ax22.scatter(final_beta_vals, stella_gamma_vals, c="black", s=mymarkersize, marker="x")

    myfile = open(gs2_pickle, 'rb')
    gs2_dict = pickle.load(myfile)  # contains 'beta', 'frequency', 'growth rate'
    myfile.close()
    gs2_beta = gs2_dict['beta']; gs2_frequency=gs2_dict['frequency']; gs2_growth_rate = gs2_dict['growth rate']
    # Sort the beta vals into ascending order
    gs2_beta = np.array(gs2_beta).flatten()
    sort_idxs = np.argsort(gs2_beta)
    gs2_beta = np.array(gs2_beta)[sort_idxs]
    gs2_frequency = np.array(gs2_frequency)[sort_idxs]
    gs2_growth_rate = np.array(gs2_growth_rate)[sort_idxs]
    ax21.plot(gs2_beta, gs2_frequency, label="GS2", lw=mylinewidth)
    ax22.plot(gs2_beta, gs2_growth_rate, label="GS2", lw=mylinewidth)

    ax21.set_ylabel(r"$\omega(a/v_{th,r})$", fontsize=myaxlabelsize)
    ax22.set_ylabel(r"$\gamma(a/v_{th,r})$", fontsize=myaxlabelsize)
    ax22.set_xlabel(r"$\beta$", fontsize=myaxlabelsize)

    ax22.legend(loc="best", fontsize=mylegendsize)
    for ax in [ax21, ax22]:
        ax.set_xlim([-0.001, 0.05])
        ax.grid(True)
        ax.tick_params(axis='both', which='major', labelsize=mylabelsize, direction="in")
        ax.tick_params(axis='both', which='minor', direction="in")

    fig2.tight_layout()
    fig2.savefig(save_name)
    plt.close(fig2)

    return

def make_ky_scan_plots(stella_longnames, gs2_longnames, ky_vals, save_name,
                         gs2_pickle=None,  plot_apar=False, plot_bpar=False, plot_format=".eps"):
    """Construct ky scans in the case where we have a set of stella sims and a
    set of GS2 sims."""


    gamma_vals = []
    freq_vals = []
    final_ky_vals = []

    stella_gamma_vals = []
    stella_freq_vals = []
    gs2_gamma_vals = []
    gs2_freq_vals = []

    for sim_idx, stella_longname in enumerate(stella_longnames):
        ## Plot of omega(t)
        try:
            make_comparison_plots([stella_longname,
                                    gs2_longnames[sim_idx]
                                    ],
                                  ["stella ky = " + str(ky_vals[sim_idx]),
                                  "GS2 ky = " + str(ky_vals[sim_idx])
                                  ],
                                  (save_name + "_ky = " + str(ky_vals[sim_idx])),
                                  sim_types=["stella",
                                            "gs2"
                                            ],
                                  plot_apar=plot_apar, plot_bpar=plot_bpar, plot_format=plot_format)
            # Currently we get gamma_stable, but don't do anything with it
            # (taking growth rate from .out is too inaccurate)
            time, freqom_final, gammaom_final, freqom, gammaom, gamma_stable = get_omega_data(stella_longname, "stella")
            freqom_final_gs2, gammaom_final_gs2 = get_gs2_omega_from_plaintext(gs2_longnames[sim_idx])
            if (np.isfinite(gammaom_final) and np.isfinite(freqom_final)
                 and np.isfinite(freqom_final_gs2) and np.isfinite(gammaom_final_gs2)
                  ):
                stella_gamma_vals.append(gammaom_final)
                stella_freq_vals.append(freqom_final)
                gs2_gamma_vals.append(gammaom_final_gs2)
                gs2_freq_vals.append(freqom_final_gs2)
                final_ky_vals.append(ky_vals[sim_idx])
        except ValueError:
            # Probably occurred because a sim has failed, and thus the .omega file,
            # or one of the other files, is bad.
            print("ValueError for sim_idx ", sim_idx)
            print("stella_longname = ", stella_longname)
            print("gs2_longname = ", gs2_longnames[sim_idx])

    stella_gamma_vals = np.array(stella_gamma_vals)
    stella_freq_vals = np.array(stella_freq_vals)
    gs2_gamma_vals = np.array(gs2_gamma_vals)
    gs2_freq_vals = np.array(gs2_freq_vals)
    final_ky_vals = np.array(final_ky_vals)
    ## Plot of omega(ky)
    fig2 = plt.figure(figsize=[10, 12])
    ax21 = fig2.add_subplot(221)
    ax22 = fig2.add_subplot(223, sharex=ax21)
    ax23 = fig2.add_subplot(222)
    ax24 = fig2.add_subplot(224, sharex=ax23)
    ax21.plot(final_ky_vals, stella_freq_vals, label="stella")
    ax22.plot(final_ky_vals, stella_gamma_vals, label="stella")
    ax21.plot(final_ky_vals, gs2_freq_vals, label="GS2")
    ax22.plot(final_ky_vals, gs2_gamma_vals, label="GS2")
    ax23.plot(final_ky_vals, (stella_freq_vals-gs2_freq_vals)/gs2_freq_vals)
    ax24.plot(final_ky_vals, (stella_gamma_vals-gs2_gamma_vals)/gs2_gamma_vals)
    ax21.scatter(final_ky_vals, stella_freq_vals, c="black", s=15., marker="x")
    ax22.scatter(final_ky_vals, stella_gamma_vals, c="black", s=15., marker="x")
    ax21.scatter(final_ky_vals, gs2_freq_vals, c="red", s=15., marker="x")
    ax22.scatter(final_ky_vals, gs2_gamma_vals, c="red", s=15., marker="x")
    if gs2_pickle is not None:
        print("Not supported")
        sys.exit()

    ax21.set_ylabel(r"$\omega$")
    ax23.set_ylabel(r"$(\Delta \omega)/\omega_{gs2}$")
    ax24.set_ylabel(r"$(\Delta \gamma)/\gamma_{gs2}$")
    ax22.set_ylabel(r"$\gamma$")
    ax22.set_xlabel(r"$k_y$")
    ax24.set_xlabel(r"$k_y$")

    for ax in [ax21, ax22]:
        ax.grid(True)
        ax.legend(loc="best")
    ax23.grid(True)
    ax24.grid(True)
    plt.tight_layout()
    fig2.savefig(save_name + "_omega_ky" + plot_format)
    plt.close(fig2)

    return

def plot_gmvus(stella_outnc_longname, which="gvpa", plot_gauss_squared=False,
                stretch_electron_vpa=False):
    """ """
    t, z, mu, vpa, gds2, gds21, gds22, bmag, gradpar, gvmus = extract_data_from_ncdf(stella_outnc_longname,
                                    "t", 'zed', "mu", "vpa", 'gds2', 'gds21', 'gds22', 'bmag', 'gradpar', 'gvmus')
    print("len(t)", "len(z), len(mu), len(vpa) = ", len(t), len(z), len(mu), len(vpa))
    #print("gvmus.shape = ", gvmus.shape)

    #sys.exit()

    if ((which == "gvpa") or (which == "both")):
        ## Code to plot g for stella
        gvmus = gvmus[-1]   # spec, mu, vpa
        fig = plt.figure(figsize=[10,10])
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex=ax1)
        counter=0
        ion_max_val = 0
        electron_max_val = 0
        for mu_idx in range(0, len(mu)):
            counter += 1
            g_ion_vpa = gvmus[0, mu_idx, :]
            g_electron_vpa = gvmus[1, mu_idx, :]
            tmp_max_val = np.max(g_ion_vpa)
            ion_max_val = max(tmp_max_val, ion_max_val)
            tmp_max_val = np.max(g_electron_vpa)
            electron_max_val = max(tmp_max_val, electron_max_val)

            ax1.plot(vpa, g_ion_vpa, label="mu={:.3f}".format(mu[mu_idx]))
            if stretch_electron_vpa:
                ax2.plot(vpa/np.sqrt(0.00028), g_electron_vpa, label="mu={:.3f}".format(mu[mu_idx]))
            else:
                ax2.plot(vpa, g_electron_vpa, label="mu={:.3f}".format(mu[mu_idx]))

            if counter == 5:
                if plot_gauss_squared:
                    maxwell_vpa_squared = (np.exp(-vpa**2*0.00028))**2
                    ax1.plot(vpa, maxwell_vpa_squared*ion_max_val, c="black", ls="--", label=r"$\exp\left(-v_\parallel^2\right)$")
                    ax2.plot(vpa, maxwell_vpa_squared*electron_max_val, c="black", ls="--", label=r"$\exp\left(-v_\parallel^2\right)$")
                for ax in [ax1, ax2]:
                    ax.grid(True)
                    ax.legend(loc="best")
                ax2.set_xlabel("vpa")
                ax1.set_ylabel(r"$g_{ion}$")
                ax2.set_ylabel(r"$g_{electron}$")
                plt.show()
                fig = plt.figure(figsize=[10,10])
                ax1 = fig.add_subplot(211)
                ax2 = fig.add_subplot(212)
                counter=0
                ion_max_val = 0
                electron_max_val = 0

        for ax in [ax1, ax2]:
            ax.grid(True)
            ax.legend(loc="best")
        ax2.set_xlabel("vpa")
        ax1.set_ylabel(r"$g_{ion}$")
        ax2.set_ylabel(r"$g_{electron}$")
        plt.show()

    if ((which == "gmu") or (which == "both")):
        ###### Plot g(mu) for different vpa
        if which == "gmu":
            # Get the final timestep of gvmus
            gvmus = gvmus[-1]   # spec, mu, vpa
        fig = plt.figure(figsize=[10,10])
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        counter=0
        #print("vpa = ", vpa)
        #sys.exit()
        # Only go over half of vpa-space (symmetry arguyment)
        for vpa_idx in range(int(len(vpa)/2), len(vpa)):
            counter += 1
            g_ion_mu = gvmus[0, :, vpa_idx]
            g_electron_mu = gvmus[1, :, vpa_idx]
            ax1.plot(mu, g_ion_mu, label="vpa={:.3f}".format(vpa[vpa_idx]))
            ax2.plot(mu, g_electron_mu, label="vpa={:.3f}".format(vpa[vpa_idx]))

            if counter == 5:
                for ax in [ax1, ax2]:
                    ax.grid(True)
                    ax.legend(loc="best")
                ax2.set_xlabel("mu")
                ax1.set_ylabel(r"$g_{ion}$")
                ax2.set_ylabel(r"$g_{electron}$")
                plt.show()
                fig = plt.figure(figsize=[10,10])
                ax1 = fig.add_subplot(211)
                ax2 = fig.add_subplot(212)
                counter=0

        for ax in [ax1, ax2]:
            ax.grid(True)
            ax.legend(loc="best")
        ax2.set_xlabel("mu")
        ax1.set_ylabel(r"$g_{ion}$")
        ax2.set_ylabel(r"$g_{electron}$")
        plt.show()

    return

def plot_gzvs(stella_outnc_longname, which="gz", plot_gauss_squared=False,
                stretch_electron_vpa=False):
    """ """
    view_ncdf_variables(stella_outnc_longname)
    t, z, vpa, gzvs = extract_data_from_ncdf(stella_outnc_longname,
                                    "t", 'zed', "vpa", 'gzvs')
    print("len(t)", "len(z), len(vpa) = ", len(t), len(z), len(vpa))
    print("gzvs.shape = ", gzvs.shape)  # t, spec, vpa, z, tube
    #sys.exit()

    if ((which == "gvpa") or (which == "both")):
        ## Code to plot g for stella
        gzvs = gzvs[-1,:,:,:,0]   # spec, vpa, z, # tube
        fig = plt.figure(figsize=[10,10])
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex=ax1)
        counter=0
        ion_max_val = 0
        electron_max_val = 0
        for z_idx in range(int(len(z)/2), len(z)):
            counter += 1
            g_ion_vpa = gzvs[0, :, z_idx]
            g_electron_vpa = gzvs[1, :, z_idx]
            tmp_max_val = np.max(g_ion_vpa)
            ion_max_val = max(tmp_max_val, ion_max_val)
            tmp_max_val = np.max(g_electron_vpa)
            electron_max_val = max(tmp_max_val, electron_max_val)

            ax1.plot(vpa, g_ion_vpa, label="z/pi={:.3f}".format(z[z_idx]/np.pi))
            ax2.plot(vpa, g_electron_vpa, label="z/pi={:.3f}".format(z[z_idx]/np.pi))

            if counter == 5:
                if plot_gauss_squared:
                    maxwell_vpa_squared = (np.exp(-vpa**2))**2
                    ax1.plot(vpa, maxwell_vpa_squared*ion_max_val, c="black", ls="--", label=r"$\exp\left(-v_\parallel^2\right)$")
                    ax2.plot(vpa, maxwell_vpa_squared*electron_max_val, c="black", ls="--", label=r"$\exp\left(-v_\parallel^2\right)$")
                for ax in [ax1, ax2]:
                    ax.grid(True)
                    ax.legend(loc="best")
                ax2.set_xlabel("vpa")
                ax1.set_ylabel(r"$g_{ion}$")
                ax2.set_ylabel(r"$g_{electron}$")
                plt.show()
                fig = plt.figure(figsize=[10,10])
                ax1 = fig.add_subplot(211)
                ax2 = fig.add_subplot(212)
                counter=0
                ion_max_val = 0
                electron_max_val = 0

        for ax in [ax1, ax2]:
            ax.grid(True)
            ax.legend(loc="best")
        ax2.set_xlabel("vpa")
        ax1.set_ylabel(r"$g_{ion}$")
        ax2.set_ylabel(r"$g_{electron}$")
        plt.show()

    if ((which == "gz") or (which == "both")):
        ###### Plot g(z) for different vpa
        if which == "gz":
            # Get the final timestep of gvmus
            gzvs = gzvs[-1,:,:,:,0]   # spec, vpa, z, # tube
            #gzvs = gzvs[:,:,:,0]   # spec, vpa, z
        fig = plt.figure(figsize=[10,10])
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        counter=0
        #print("vpa = ", vpa)
        #sys.exit()
        # Only go over half of vpa-space (symmetry argument)
        for vpa_idx in range(int(len(vpa)/2), len(vpa)):
            counter += 1
            g_ion_z = gzvs[0, vpa_idx, :]
            g_electron_z = gzvs[1, vpa_idx, :]
            ax1.plot(z/np.pi, g_ion_z, label="vpa={:.3f}".format(vpa[vpa_idx]))
            ax2.plot(z/np.pi, g_electron_z, label="vpa={:.3f}".format(vpa[vpa_idx]))

            if counter == 5:
                for ax in [ax1, ax2]:
                    ax.grid(True)
                    ax.legend(loc="best")
                ax2.set_xlabel(r"$z/\pi$")
                ax1.set_ylabel(r"$g_{ion}$")
                ax2.set_ylabel(r"$g_{electron}$")
                plt.show()
                fig = plt.figure(figsize=[10,10])
                ax1 = fig.add_subplot(211)
                ax2 = fig.add_subplot(212)
                counter=0

        for ax in [ax1, ax2]:
            ax.grid(True)
            ax.legend(loc="best")
        ax2.set_xlabel(r"$z/\pi$")
        ax1.set_ylabel(r"$g_{ion}$")
        ax2.set_ylabel(r"$g_{electron}$")
        plt.show()

    return
