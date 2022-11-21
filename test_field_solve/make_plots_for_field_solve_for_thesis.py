"""" """

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
from scipy.special import iv
import make_sims_for_thesis as make_sims
import os

gs2_sim_fapar0_fbpar0 = "sims/gs2_slab_fapar0_fbpar0/input.out.nc"
gs2_sim_fapar0_fbpar0_lr = "sims/gs2_slab_fapar0_fbpar0/input_ngauss6_negrid8.out.nc"
gs2_sim_fapar0_fbpar0_vlr = "sims/gs2_slab_fapar0_fbpar0/input_ngauss3_negrid4.out.nc"
gs2_sim_fapar1_fbpar1 = "sims/gs2_slab_fapar1_fbpar1/input.out.nc"
gs2_sim_fapar1_fbpar1_lr = "sims/gs2_slab_fapar1_fbpar1/input_ngauss6_negrid8.out.nc"
gs2_sim_fapar1_fbpar1_vlr = "sims/gs2_slab_fapar1_fbpar1/input_ngauss3_negrid4.out.nc"
gs2_sim_fapar0_fbpar1 = "sims/gs2_slab_fapar0_fbpar1/input.out.nc"
gs2_sim_fapar1_fbpar0 = "sims/gs2_slab_fapar1_fbpar0/input.out.nc"
stella_sim_fapar0_fbpar0 = "sims/stella_slab_fapar0_fbpar0/input.out.nc"
stella_sim_fapar0_fbpar0_lr = "sims/stella_slab_fapar0_fbpar0/input_nvgrid12_nmu6.out.nc"
stella_sim_fapar0_fbpar0_vlr = "sims/stella_slab_fapar0_fbpar0/input_nvgrid6_nmu3.out.nc"
stella_sim_fapar1_fbpar1_h= "sims/stella_slab_fapar1_fbpar1_h/input.out.nc"
stella_sim_fapar1_fbpar1= "sims/stella_slab_fapar1_fbpar1/input.out.nc"
stella_sim_fapar1_fbpar1_lr= "sims/stella_slab_fapar1_fbpar1/input_nvgrid12_nmu6.out.nc"
stella_sim_fapar1_fbpar1_vlr= "sims/stella_slab_fapar1_fbpar1/input_nvgrid6_nmu3.out.nc"
stella_sim_fapar0_fbpar1= "sims/stella_slab_fapar0_fbpar1/input.out.nc"
stella_sim_fapar1_fbpar0= "sims/stella_slab_fapar1_fbpar0/input.out.nc"

def compare_field_for_thesis(field_name, outnc_longnames, sim_types, sim_labels,
                            # normalise=False,
                            save_name="a.eps",
                            ax1_yticks=None,
                            ax1_yticklabels=None,
                            ax2_yticks=None,
                            ax2_yticklabels=None):
    """ """

    ylabel_fontsize = 30
    xlabel_fontsize = 30
    legend_fontsize = 18
    linestyles=((0,(1,0)), (0, (2,2)), (0,(5,1,1,1)), (0,(1,1)),
                (0, (3,2,2,3)), (0, (1,0)) )
    # linewidths = [6, 4, 3, 3, 3, 2]

    yticklabel_fontsize = 25
    xticklabel_fontsize = yticklabel_fontsize
    my_xticklength = 5
    my_xtickwidth = 1
    my_yticklength = 5
    my_ytickwidth = 1

    default_cols = plt.get_cmap("tab10")

    left = 0.15
    right = 0.985
    top = 0.985
    bottom = 0.1
    vspace=0.05
    height = (top - bottom - vspace)/2
    width = right - left
    my_linewidth = 3

    if (field_name != "phi") and (field_name != "apar") and (field_name != "bpar"):
        print("field_name not recognised!")
        sys.exit()
    ### Keys for gs2, stella
    field_key_stella = field_name + "_vs_t"
    field_key_gs2 = field_name + "_t"
    ### Initialise lists to store the data
    z_list = [] ; field_real_list = [] ; field_imag_list = []
    sim_labels_list = []
    ### Try to get the data
    for sim_idx, outnc_longname in enumerate(outnc_longnames):
        try:
            if sim_types[sim_idx] == "stella":
                (z, field_vs_t) = extract_data_from_ncdf_with_xarray(outnc_longname,
                        'zed', field_key_stella)
                field_t0 = np.array(field_vs_t[0,0,:,0,0,:])
                # shape is now (zed, ri)
            elif sim_types[sim_idx] == "gs2":
                (z, field_vs_t) = extract_data_from_ncdf_with_xarray(outnc_longname,
                        'theta', field_key_gs2)
                field_t0 = np.array(field_vs_t[0,0,0,:,:])
                # shape is now (zed, ri)
            else:
                print("sim_type not recognised!")
                sys.exit()
            z_list.append(z)
            field_real_list.append(np.array(field_t0[:,0]))
            field_imag_list.append(np.array(field_t0[:,1]))
            sim_labels_list.append(sim_labels[sim_idx])
        except KeyError:
            print("field not found. outnc_longname, field = ", outnc_longname, field_name)

    fig = plt.figure(figsize=(12,12))
    ax1 = fig.add_axes((left, bottom + height + vspace, width, height))
    ax2 = fig.add_axes((left, bottom, width, height))

    z_orig = np.array(z_list[0])
    field_real = np.array(field_real_list[0])
    for sim_idx, sim_label in enumerate(sim_labels_list):
        z = np.array(z_list[sim_idx])
        field_real = np.array(field_real_list[sim_idx])
        field_imag = np.array(field_imag_list[sim_idx])
        # if normalise:
        #     norm_fac_real = max(np.max(abs(field_real)), 1E-12)
        #     field_real = field_real/norm_fac_real
        #     norm_fac_imag = max(np.max(abs(field_imag)), 1E-12)
        #     field_imag = field_imag/norm_fac_imag
        ax1.plot(z/np.pi, field_real, ls=linestyles[sim_idx], lw=my_linewidth, label=sim_label, c=default_cols(sim_idx))

        if len(z_orig) == len(z):
            if np.max(abs(z_orig - z)) < 1E-4:
                if sim_idx != 0:    # Don't plot the first, as this is what we're compaing to
                    ax2.plot(z/np.pi, (abs((field_real - field_real_list[0])/field_real_list[0])*100),
                                        c=default_cols(sim_idx), ls=linestyles[sim_idx], lw=my_linewidth)
            else:
                print("z vals don't match")
                print("z_orig - z = ", z_orig - z)
                print("np.max(abs(z_orig - z)) = ", np.max(abs(z_orig - z)))
                print("z_orig = ", np.array(z_orig))
                print("z = ", np.array(z))
        else:
            print("z lengths not equal")

    ax2.set_yscale("log")

    if field_name == "phi":
        field_ylabel_ax1 = r"$\tilde{\phi}_k$"
        field_ylabel_ax2 = r"$\tilde{\phi}_k$ error (%)"
    if field_name == "apar":
        field_ylabel_ax1 = r"$\tilde{A}_{\parallel, k}$"
        field_ylabel_ax2 = r"$\tilde{A}_{\parallel, k}$ error (%)"
    if field_name == "bpar":
        field_ylabel_ax1 = r"$\tilde{B}_{\parallel, k}$"
        field_ylabel_ax2 = r"$\tilde{B}_{\parallel, k}$ error (%)"

    ax1.set_xlim((-1, 1))
    ax2.set_xlim((-1, 1))
    if ax1_yticks is not None:
        ax1.set_yticks(ax1_yticks)
        ax1.set_yticklabels(ax1_yticklabels, fontsize=yticklabel_fontsize)
    if ax2_yticks is not None:
        ax2.set_yticks(ax2_yticks)
        ax2.set_yticklabels(ax2_yticklabels, fontsize=yticklabel_fontsize)
    ax1.tick_params("x", length=my_xticklength, width=my_xtickwidth, direction="out")
    ax2.tick_params("x", length=my_xticklength, width=my_xtickwidth, direction="out")
    ax1.tick_params("y", length=my_yticklength, width=my_ytickwidth, direction="out")
    ax1.set_xticks([-1, 0, 1])
    ax1.set_xticklabels([])
    ax2.set_xticks([-1, 0, 1])
    ax2.set_xticklabels(["-1", "0", "1"], fontsize=xticklabel_fontsize)
    ax1.legend(loc="best", fontsize=legend_fontsize)
    ax1.set_ylabel(field_ylabel_ax1, fontsize=ylabel_fontsize)
    ax2.set_ylabel(field_ylabel_ax2, fontsize=ylabel_fontsize, labelpad=20)
    ax2.set_xlabel(r"$z/\pi$", fontsize=xlabel_fontsize)

    plt.savefig(save_name)
    plt.close()
    return

def compare_phi_for_thesis(outnc_longnames, sim_types, sim_labels, **kwargs):
    """ """
    compare_field_for_thesis("phi", outnc_longnames, sim_types, sim_labels, **kwargs)
    return

def compare_apar_for_thesis(outnc_longnames, sim_types, sim_labels, **kwargs):
    """ """
    compare_field_for_thesis("apar", outnc_longnames, sim_types, sim_labels, **kwargs)
    return

def compare_bpar_for_thesis(outnc_longnames, sim_types, sim_labels, **kwargs):
    """ """
    compare_field_for_thesis("bpar", outnc_longnames, sim_types, sim_labels, **kwargs)
    return

def calculate_fields_zjzeroexp(kperp, beta, dist, m_ion=1, m_electron=2.8E-4, B=1):
    """ """
    b_ion = 0.5*kperp*kperp*m_ion/(B*B)
    b_electron = 0.5*kperp*kperp*m_electron/(B*B)
    gamzero_ion = np.exp(-b_ion)*iv(0,b_ion)
    gamzero_electron = np.exp(-b_electron)*iv(0,b_electron)
    gamone_ion = np.exp(-b_ion)*(iv(0,b_ion) - iv(1,b_ion) )
    gamone_electron = np.exp(-b_electron)*(iv(0,b_electron) - iv(1,b_electron) )
    antot1 = (gamzero_ion + gamzero_electron)
    antot3 = (-0.5*beta/B) * (gamone_ion - gamone_electron )
    # print("gamzero_ion = ", gamzero_ion)
    # print("gamzero_electron = ", gamzero_electron)
    if dist == "h":
        gamtot_h = 2
        phi = antot1/gamtot_h
        apar = 0
        bpar = antot3
    else:
        gamtot = 2 - (gamzero_ion + gamzero_electron)
        gamtot13 =  (gamone_electron - gamone_ion)/B
        gamtot31 = beta*(gamone_ion - gamone_electron)/(2*B)
        gamtot33 = 1 + beta*(gamone_ion + gamone_electron)/B

        phi = (antot1 - (gamtot13/gamtot33)*antot3 ) / (gamtot - gamtot13*gamtot31/gamtot33 )
        apar = 0
        bpar = (antot3 - (gamtot31/gamtot)*antot1) / (gamtot33 - gamtot13*gamtot31/gamtot)
    return phi, apar, bpar

def calculate_fields_zvpajzeroexp(kperp, beta, dist, m_ion=1, m_electron=2.8E-4, B=1):
    """ """
    b_ion = 0.5*kperp*kperp*m_ion/(B*B)
    b_electron = 0.5*kperp*kperp*m_electron/(B*B)
    gamzero_ion = np.exp(-b_ion)*iv(0,b_ion)
    gamzero_electron = np.exp(-b_electron)*iv(0,b_electron)
    gamone_ion = np.exp(-b_ion)*(iv(0,b_ion) - iv(1,b_ion) )
    gamone_electron = np.exp(-b_electron)*(iv(0,b_electron) - iv(1,b_electron) )
    antot2 = 0.5*beta*(gamzero_ion/np.sqrt(m_ion) + gamzero_electron/np.sqrt(m_electron))
    # print("np.exp(-b_ion), iv(0,b_ion) = ", np.exp(-b_ion), iv(0,b_ion))
    # print("gamzero_ion, gamzero_electron, gamone_ion, gamone_electron = ", gamzero_ion, gamzero_electron, gamone_ion, gamone_electron)
    if dist == "h":
        apar_denom_h = kperp*kperp
        phi = 0
        apar = antot2/apar_denom_h
        bpar = 0
    else:
        apar_denom = kperp*kperp + beta*(gamzero_ion/m_ion + gamzero_electron/m_electron)
        phi = 0
        apar = antot2/apar_denom
        bpar = 0
    return phi, apar, bpar

def get_fields_from_outnc_files(folder_longname, kperp=False):
    """ """
    ## If the pickle exists, use the pickle
    pickle_longname = "sims/" + folder_longname + "/field_results.pickle"
    print("folder_longname = ", folder_longname)
    if os.path.isfile(pickle_longname):
        myfile = open(pickle_longname, "rb")
        [vals_1, vals_2, phi_vals, apar_vals, bpar_vals] = pickle.load(myfile)
        myfile.close()
    else:
        outnc_longnames = glob.glob("sims/" + folder_longname + "/*.out.nc")
        print("outnc_longnames = ", outnc_longnames)
        vals_1 = []
        vals_2 = []
        phi_vals = []
        apar_vals = []
        bpar_vals = []
        for outnc_longname in outnc_longnames:
            sim_longname = re.split(".out.nc", outnc_longname)[0]
            vals_str_list = re.split("_", sim_longname)
            if kperp:
                vals_1.append(float(vals_str_list[-3]))
            else:
                vals_1.append(float(vals_str_list[-2]))
                vals_2.append(int(vals_str_list[-1]))
            field_key_stella = "phi_vs_t"
            (phi_vs_t, apar_vs_t, bpar_vs_t) = extract_data_from_ncdf_with_xarray(outnc_longname,
                            "phi_vs_t", "apar_vs_t", "bpar_vs_t")
            phi_vals.append(phi_vs_t[0,0,0,0,0,0])
            apar_vals.append(apar_vs_t[0,0,0,0,0,0])
            bpar_vals.append(bpar_vs_t[0,0,0,0,0,0])
        vals_1 = np.array(vals_1)
        vals_2 = np.array(vals_2)
        phi_vals = np.array(phi_vals)
        apar_vals = np.array(apar_vals)
        bpar_vals = np.array(bpar_vals)

        myfile = open(pickle_longname, "wb")
        pickle.dump([vals_1, vals_2, phi_vals, apar_vals, bpar_vals], myfile)
        myfile.close()

    return vals_1, vals_2, phi_vals, apar_vals, bpar_vals

def vspace_res_test_field_solve_for_thesis():
    """Plot % error in phi, bpar and abs(apar) in various cases:
    (1) Fixed vperpmax, nmu. plotting (error in field) vs nvpa for several vpamax
    (2) Fixed vpamax, nvpa, plotting (error in field) vs nvmu for several vperpmax
    For each of these, plot the result for field solve in h and field solve in h.
    """

    marker_size = 120
    legend_fontsize = 13
    marker_list = ["s", "o", "P", "X", "v", "^", "<", ">", "1", "2", "3"]
    lw_list = [4, 3, 3, 2]
    ls_list = ["-", "--", "-.", (0, (4,1,2,1))]
    title_fontsize = 40
    my_xticklength = 7
    my_xtickwidth = 2
    my_xminorticklength = 4
    my_xminortickwidth = 1
    my_yticklength = 7
    my_ytickwidth = 2
    my_yminorticklength = 4
    my_yminortickwidth = 1

    xticklabel_fontsize = 20
    yticklabel_fontsize = 20
    ylabel_fontsize = 36
    xlabel_fontsize = 36
    bracket_fontsize = 60

    top = 0.94
    left = 0.14
    right = 0.98
    bottom = 0.07
    vspace = 0.02
    hspace = 0.08

    def init_plot():
        """ """

        height = (top - bottom - 2*vspace)/3
        row3_bottom = bottom
        row2_bottom = bottom + height + vspace
        row1_bottom = row2_bottom + height + vspace
        width = (right - left - hspace)/2
        col1_left = left
        col2_left = left + width + hspace

        fig = plt.figure(figsize=(12,12))
        ax1 = fig.add_axes((col1_left, row1_bottom, width, height))
        ax2 = fig.add_axes((col2_left, row1_bottom, width, height))
        ax3 = fig.add_axes((col1_left, row2_bottom, width, height))
        ax4 = fig.add_axes((col2_left, row2_bottom, width, height))
        ax5 = fig.add_axes((col1_left, row3_bottom, width, height))
        ax6 = fig.add_axes((col2_left, row3_bottom, width, height))

        return fig, ax1, ax2, ax3, ax4, ax5, ax6

    def make_plot(fig, ax1, ax2, ax3, ax4, ax5, ax6, which):
        """ """
        if which=="vpa":
            unique_vpamax = sorted(set(vpamax_vals_h))
            vpamax_array = np.array(vpamax_vals_h)
            nvpa_array = np.array(nvpa_vals_h)
            # Calculate the magnitude of % error for phi, bpar
            # and magnitude of error for apar (not %, since analytic apar = 0)
            phi_array = np.abs(100*((phi_vals_h) - analytic_phi_h)/analytic_phi_h)
            apar_array = np.abs((apar_vals_h))
            bpar_array = np.abs(100*((bpar_vals_h) - analytic_bpar_h)/analytic_bpar_h)

            for vpamax_idx, unique_vpamax_val in enumerate(unique_vpamax):
                idxs = np.argwhere(np.abs(vpamax_array-unique_vpamax_val)<1E-4)
                ax1.scatter(nvpa_array[idxs], phi_array[idxs], marker=marker_list[vpamax_idx], s=marker_size)
                ax3.scatter(nvpa_array[idxs], apar_array[idxs], marker=marker_list[vpamax_idx],
                            s=marker_size, label=r"$\tilde{v}_{\parallel, \rm max}$="+str(unique_vpamax_val))
                ax5.scatter(nvpa_array[idxs], bpar_array[idxs], marker=marker_list[vpamax_idx], s=marker_size)

            unique_vpamax = sorted(set(vpamax_vals_gbar))
            vpamax_array = np.array(vpamax_vals_gbar)
            nvpa_array = np.array(nvpa_vals_gbar)
            # Calculate the magnitude of % error for phi, bpar
            # and magnitude of error for apar (not %, since analytic apar = 0)
            phi_array = np.abs(100*((phi_vals_gbar) - analytic_phi_gbar)/analytic_phi_gbar)
            apar_array = np.abs((apar_vals_gbar))
            bpar_array = np.abs(100*((bpar_vals_gbar) - analytic_bpar_gbar)/analytic_bpar_gbar)
            for vpamax_idx, unique_vpamax_val in enumerate(unique_vpamax):
                idxs = np.argwhere(np.abs(vpamax_array-unique_vpamax_val)<1E-4)
                ax2.scatter(nvpa_array[idxs], phi_array[idxs], marker=marker_list[vpamax_idx], s=marker_size)
                ax4.scatter(nvpa_array[idxs], apar_array[idxs], marker=marker_list[vpamax_idx],
                            s=marker_size, label=r"$\tilde{v}_{\parallel, \rm max}$="+str(unique_vpamax_val))
                ax6.scatter(nvpa_array[idxs], bpar_array[idxs], marker=marker_list[vpamax_idx], s=marker_size)

        elif which=="vperp":
            unique_vperpmax = sorted(set(vperpmax_vals_h))
            vperpmax_array = np.array(vperpmax_vals_h)
            # Calculate the magnitude of % error for phi, bpar
            # and magnitude of error for apar (not %, since analytic apar = 0)
            nmu_array = np.array(nmu_vals_h)
            phi_array = np.abs(100*((phi_vals_h) - analytic_phi_h)/analytic_phi_h)
            apar_array = np.abs((apar_vals_h))
            bpar_array = np.abs(100*((bpar_vals_h) - analytic_bpar_h)/analytic_bpar_h)

            for vperpmax_idx, unique_vperpmax_val in enumerate(unique_vperpmax):
                idxs = np.argwhere(np.abs(vperpmax_array-unique_vperpmax_val)<1E-4)
                ax1.scatter(nmu_array[idxs], phi_array[idxs], marker=marker_list[vperpmax_idx], s=marker_size)
                ax3.scatter(nmu_array[idxs], apar_array[idxs], marker=marker_list[vperpmax_idx],
                            s=marker_size, label=r"$\tilde{v}_{\perp, \rm max}$="+str(unique_vperpmax_val))
                ax5.scatter(nmu_array[idxs], bpar_array[idxs], marker=marker_list[vperpmax_idx], s=marker_size)
            ## Repeat for gbar dist.
            unique_vperpmax = sorted(set(vperpmax_vals_gbar))
            vperpmax_array = np.array(vperpmax_vals_gbar)
            # Calculate the magnitude of % error for phi, bpar
            # and magnitude of error for apar (not %, since analytic apar = 0)
            nmu_array = np.array(nmu_vals_gbar)
            phi_array = np.abs(100*(np.array(phi_vals_gbar) - analytic_phi_gbar)/analytic_phi_gbar)
            apar_array = np.abs(np.array(apar_vals_gbar))
            bpar_array = np.abs(100*(np.array(bpar_vals_gbar) - analytic_bpar_gbar)/analytic_bpar_gbar)

            for vperpmax_idx, unique_vperpmax_val in enumerate(unique_vperpmax):
                idxs = np.argwhere(np.abs(vperpmax_array-unique_vperpmax_val)<1E-4)
                ax2.scatter(nmu_array[idxs], phi_array[idxs], marker=marker_list[vperpmax_idx], s=marker_size)
                ax4.scatter(nmu_array[idxs], apar_array[idxs], marker=marker_list[vperpmax_idx],
                            s=marker_size, label=r"$\tilde{v}_{\perp, \rm max}$="+str(unique_vperpmax_val))
                ax6.scatter(nmu_array[idxs], bpar_array[idxs], marker=marker_list[vperpmax_idx], s=marker_size)

        return

    def finish_plot(fig, ax1, ax2, ax3, ax4, ax5, ax6, which):
        """ """
        ########## Beginning stuff
        ax1.set_title(r"$\tilde{h}_{k}$", fontsize=title_fontsize)
        ax2.set_title(r"$\tilde{\bar{g}}_{k}$", fontsize=title_fontsize)
        for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
            #ax.set_xlim((0, 150))
            ax.set_xscale("log")
            ax.set_yscale("log")
        ax1.set_ylabel(r"$\vert\Delta \tilde{\phi}_k\vert (\%)$ ", fontsize=ylabel_fontsize)
        ax3.set_ylabel(r"$\vert\Delta \tilde{A}_{\parallel k}\vert$ ", fontsize=ylabel_fontsize)
        ax5.set_ylabel(r"$\vert\Delta \tilde{B}_{\parallel k}\vert (\%)$ ", fontsize=ylabel_fontsize)
        ax3.legend(loc="lower right", fontsize=legend_fontsize, ncol=2)

        for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
            ax.set_xticks([10, 100])
            ax.tick_params("x", which="major", length=my_xticklength, width=my_xtickwidth, direction="in")
            ax.tick_params("x", which="minor", length=my_xminorticklength, width=my_xminortickwidth, direction="in")
            ax.tick_params("y", which="major", length=my_yticklength, width=my_ytickwidth, direction="in")
            ax.tick_params("y", which="minor", length=my_yminorticklength, width=my_yminortickwidth, direction="in")
        for ax in [ax1, ax2, ax3, ax4]:
            ax.set_xticklabels([])
        for ax in [ax5, ax6]:
            ax.set_xticklabels([r"$10^1$", r"$10^2$"], fontsize=xticklabel_fontsize)
        ########## plot-specific stuff

        if which=="vpa":
            for ax in [ax5, ax6]:
                ax.set_xlabel(r"$n_{v \parallel}$", fontsize=xlabel_fontsize)

            ax1.set_ylim([2e-13, 1e-11])
            ax1.set_yticks([1e-12, 1e-11])
            ax1.set_yticklabels([r"$10^{-12}$",
                r"$10^{-11}$"], fontsize=yticklabel_fontsize)
            # ax1.set_yticks([1e-10,1e-6, 1e-2, 100], minor=True)
            # ax1.set_yticklabels([], minor=True)

            ax2.set_yticks([1e-8, 1e-4, 1])
            ax2.set_yticklabels([r"$10^{-8}$",
                r"$10^{-4}$", r"$1$"], fontsize=yticklabel_fontsize)
            ax2.set_yticks([1e-10,1e-6, 1e-2,], minor=True)
            ax2.set_yticklabels([], minor=True)

            ax3.set_ylim([1E-19, 1.4E-15])
            ax3.set_yticks([1e-19, 1e-18, 1e-17, 1e-16, 1e-15])
            ax3.set_yticklabels([
                r"$10^{-19}$", r"$10^{-18}$", r"$10^{-17}$", r"$10^{-16}$", r"$10^{-15}$"],
                fontsize=yticklabel_fontsize)
            # ax3.set_yticks([1e-10,1e-6, 1e-2, 100], minor=True)
            # ax3.set_yticklabels([], minor=True)

            ax4.set_ylim([1E-21, 1E-18])
            ax4.set_yticks([1e-21, 1e-20, 1e-19, 1e-18])
            ax4.set_yticklabels([r"$10^{-21}$",
                r"$10^{-20}$", r"$10^{-19}$", r"$10^{-18}$"],
                fontsize=yticklabel_fontsize)
            # ax4.set_yticks([1e-10,1e-6, 1e-2, 100], minor=True)
            # ax4.set_yticklabels([], minor=True)

            ax5.set_ylim([4E-10, 160])
            ax5.set_yticks([1e-8, 1e-4, 1])
            ax5.set_yticklabels([r"$10^{-8}$",
                r"$10^{-4}$", r"$1$"], fontsize=yticklabel_fontsize)
            ax5.set_yticks([1e-10, 1e-6, 1e-2, 100], minor=True)
            ax5.set_yticklabels([], minor=True)

            ax6.set_ylim([4E-10, 160])
            ax6.set_yticks([1e-8, 1e-4, 1])
            ax6.set_yticklabels([r"$10^{-8}$",
                r"$10^{-4}$", r"$1$", ], fontsize=yticklabel_fontsize)
            ax6.set_yticks([1e-6,1e-2, 100], minor=True)
            ax6.set_yticklabels([], minor=True)

            plt.savefig("vpa_res_field_solve.eps")

        if which=="vperp":
            for ax in [ax5, ax6]:
                ax.set_xlabel(r"$n_{\mu}$", fontsize=xlabel_fontsize)

            ax1.set_yticks([1e-12,1e-8, 1e-4, 1])
            ax1.set_yticklabels([r"$10^{-12}$",r"$10^{-8}$",
                r"$10^{-4}$", r"$1$"], fontsize=yticklabel_fontsize)
            ax1.set_yticks([1e-10,1e-6, 1e-2, 100], minor=True)
            ax1.set_yticklabels([], minor=True)

            ax2.set_yticks([1e-8, 1e-4, 1])
            ax2.set_yticklabels([r"$10^{-8}$",
                r"$10^{-4}$", r"$1$"], fontsize=yticklabel_fontsize)
            ax2.set_yticks([1e-10,1e-6, 1e-2, 100], minor=True)
            ax2.set_yticklabels([], minor=True)

            ax3.set_ylim([1E-19, 3E-15])
            ax3.set_yticks([1e-19, 1e-18, 1e-17, 1e-16, 1e-15])
            ax3.set_yticklabels([r"$10^{-19}$",
                r"$10^{-18}$", r"$10^{-17}$", r"$10^{-16}$", r"$10^{-15}$"],
                fontsize=yticklabel_fontsize)
            # ax3.set_yticks([1e-10,1e-6, 1e-2, 100], minor=True)
            # ax3.set_yticklabels([], minor=True)

            ax4.set_ylim([1E-21, 1E-18])
            ax4.set_yticks([1e-21, 1e-20, 1e-19, 1e-18])
            ax4.set_yticklabels([r"$10^{-21}$",
                r"$10^{-20}$", r"$10^{-19}$", r"$10^{-18}$"],
                fontsize=yticklabel_fontsize)
            # ax4.set_yticks([1e-10,1e-6, 1e-2, 100], minor=True)
            # ax4.set_yticklabels([], minor=True)

            ax5.set_yticks([1e-8, 1e-4, 1])
            ax5.set_yticklabels([r"$10^{-8}$",
                r"$10^{-4}$", r"$1$"], fontsize=yticklabel_fontsize)
            ax5.set_yticks([1e-10, 1e-6, 1e-2, 100], minor=True)
            ax5.set_yticklabels([], minor=True)

            ax6.set_yticks([1e-8, 1e-4, 1])
            ax6.set_yticklabels([r"$10^{-8}$",
                r"$10^{-4}$", r"$1$", ], fontsize=yticklabel_fontsize)
            ax6.set_yticks([1e-10, 1e-6,1e-2, 100], minor=True)
            ax6.set_yticklabels([], minor=True)

            plt.savefig("vperp_res_field_solve.eps")

        plt.close()

        return

    analytic_phi_h, analytic_apar_h, analytic_bpar_h = calculate_fields_zjzeroexp(1, 1, "h")
    analytic_phi_gbar, analytic_apar_gbar, analytic_bpar_gbar = calculate_fields_zjzeroexp(1, 1, "gbar")
    print("h")
    ### vpa parameter tests
    vpamax_vals_h, nvpa_vals_h, phi_vals_h, apar_vals_h, bpar_vals_h = get_fields_from_outnc_files(make_sims.phi_bpar_h_vpa_scan_folder)
    vpamax_vals_gbar, nvpa_vals_gbar, phi_vals_gbar, apar_vals_gbar, bpar_vals_gbar = get_fields_from_outnc_files(make_sims.phi_bpar_gbar_vpa_scan_folder)
    fig, ax1, ax2, ax3, ax4, ax5, ax6 = init_plot()
    make_plot(fig, ax1, ax2, ax3, ax4, ax5, ax6, "vpa")
    finish_plot(fig, ax1, ax2, ax3, ax4, ax5, ax6, "vpa")

    ### vperp parameter tests
    vperpmax_vals_h, nmu_vals_h, phi_vals_h, apar_vals_h, bpar_vals_h = get_fields_from_outnc_files(make_sims.phi_bpar_h_vperp_scan_folder)
    vperpmax_vals_gbar, nmu_vals_gbar, phi_vals_gbar, apar_vals_gbar, bpar_vals_gbar = get_fields_from_outnc_files(make_sims.phi_bpar_gbar_vperp_scan_folder)
    fig, ax1, ax2, ax3, ax4, ax5, ax6 = init_plot()
    make_plot(fig, ax1, ax2, ax3, ax4, ax5, ax6, "vperp")
    finish_plot(fig, ax1, ax2, ax3, ax4, ax5, ax6, "vperp")

    return

def test_field_solve_apar_for_thesis():
    """Plot % error in apar in various cases:
    (1) Fixed vperpmax, nmu. plotting (error in field) vs nvpa for several vpamax
    (2) Fixed vpamax, nvpa, plotting (error in field) vs nvmu for several vperpmax
    For each of these, plot the result for field solve in h and field solve in h.
    (3) Fixed vspace-res, varying ky
    """

    def make_analytic_apar_plot():
        """ """
        top = 0.97
        left = 0.17
        right = 0.98
        bottom = 0.15
        vspace = 0.05
        hspace = 0.04
        height = (top - bottom)
        width = (right - left )
        col1_left = left
        col2_left = left + width + hspace

        my_xticklength = 7
        my_xtickwidth = 2
        my_xminorticklength = 4
        my_xminortickwidth = 1
        my_yticklength = 7
        my_ytickwidth = 2
        my_yminorticklength = 4
        my_yminortickwidth = 1

        label_fontsize = 40
        legend_fontsize = 50
        marker_size = 30
        xticklabel_fontsize = 30
        yticklabel_fontsize = 30

        fig = plt.figure(figsize=(12,8))
        ax1 = fig.add_axes((col1_left, bottom, width, height)) # vpa-res, h
        #ax2 = fig.add_axes((col2_left, bottom, width, height)) # vpa-res, gbar

        ax1.plot(kperp_vals, analytic_apar_h_kperp, lw=3, marker="o", mfc="none", markersize=marker_size, label=r"$\tilde{h}_k$")
        ax1.plot(kperp_vals, analytic_apar_gbar_kperp, lw=3, marker="s", mfc="none", markersize=marker_size, label=r"$\tilde{\bar{g}}_k$")
        ax1.set_xscale("log")
        ax1.legend(loc="best", fontsize=legend_fontsize)
        #ax2.set_xscale("log")
        ax1.set_yscale("log")
        ax1.set_xlabel(r"$\tilde{k}_\perp$", fontsize=label_fontsize)
        ax1.set_ylabel(r"$\tilde{A}_{\parallel k}$ (analytic)", fontsize=label_fontsize)

        ax1.tick_params("x", which="major", length=my_xticklength, width=my_xtickwidth, direction="in")
        ax1.tick_params("x", which="minor", length=my_xminorticklength, width=my_xminortickwidth, direction="in")
        ax1.tick_params("y", which="major", length=my_yticklength, width=my_ytickwidth, direction="in")
        ax1.tick_params("y", which="minor", length=my_yminorticklength, width=my_yminortickwidth, direction="in")
        ax1.set_xticks([1e-4, 1e-3, 1e-2, 1e-1, 1, 10])
        ax1.set_xticklabels([r"$10^{-4}$", r"$10^{-3}$", r"$10^{-2}$", r"$10^{-1}$", r"$1$", r"$10$"], fontsize=xticklabel_fontsize)
        ax1.set_yticks([1e-2, 1, 1e2, 1e4, 1e6, 1e8, 1e10])
        ax1.set_yticklabels([r"$10^{-2}$", r"$1$", r"$10^{2}$", r"$10^{4}$", r"$10^6$", r"$10^8$", r"$10^{10}$"], fontsize=yticklabel_fontsize)
        ax1.set_yticks([1e-1, 10, 1e3, 1e5, 1e7, 1e9], minor=True)
        ax1.set_yticklabels([], minor=True)
        #ax2.set_ylabel(r"analytic $\tilde{A}_{\parallel k}$")
        plt.savefig("analytic_apar_for_test.eps")
        plt.close()
        return

    def make_plot():
        """ """
        marker_size = 120

        legend_fontsize = 13

        my_xticklength = 7
        my_xtickwidth = 2
        my_xminorticklength = 4
        my_xminortickwidth = 1
        my_yticklength = 7
        my_ytickwidth = 2
        my_yminorticklength = 4
        my_yminortickwidth = 1

        xticklabel_fontsize = 20
        yticklabel_fontsize = 20
        ylabel_fontsize = 36
        xlabel_fontsize = 36
        bracket_fontsize = 60

        marker_list = ["s", "o", "P", "X", "v", "^", "<", ">", "1", "2", "3"]
        lw_list = [4, 3, 3, 2]
        ls_list = ["-", "--", "-.", (0, (4,1,2,1))]

        top = 0.98
        left = 0.12
        right = 0.98
        bottom = 0.13
        vspace = 0.05
        hspace = 0.06
        height = (top - bottom - vspace)/2
        row2_bottom = bottom
        row1_bottom = row2_bottom + height + vspace
        width = (right - left - 2*hspace)/3
        col1_left = left
        col2_left = left + width + hspace
        col3_left = col2_left + width + hspace

        fig = plt.figure(figsize=(14,8))
        ax1 = fig.add_axes((col1_left, row1_bottom, width, height)) # vpa-res, h
        ax2 = fig.add_axes((col1_left, row2_bottom, width, height)) # vpa-res, gbar
        ax3 = fig.add_axes((col2_left, row1_bottom, width, height)) # vperp-res, h
        ax4 = fig.add_axes((col2_left, row2_bottom, width, height)) # vperp-res, gbar
        ax5 = fig.add_axes((col3_left, row1_bottom, width, height)) # kperp, h
        ax6 = fig.add_axes((col3_left, row2_bottom, width, height)) # kperp, gbar

        ## Plot vpa res
        unique_vpamax_h = sorted(set(vpamax_vals_h))
        unique_vpamax_gbar = sorted(set(vpamax_vals_gbar))
        vpamax_array_h = np.array(vpamax_vals_h)
        nvpa_array_h = np.array(nvpa_vals_h)
        vpamax_array_gbar = np.array(vpamax_vals_gbar)
        nvpa_array_gbar = np.array(nvpa_vals_gbar)

        unique_vperpmax_h = sorted(set(vperpmax_vals_h))
        vperpmax_array_h = np.array(vperpmax_vals_h)
        nmu_array_h = np.array(nmu_vals_h)
        unique_vperpmax_gbar = sorted(set(vperpmax_vals_gbar))
        vperpmax_array_gbar = np.array(vperpmax_vals_gbar)
        nmu_array_gbar = np.array(nmu_vals_gbar)

        # Calculate the magnitude of % error for phi, bpar
        # and magnitude of error for apar (not %, since analytic apar = 0)
        phi_array_vpa = np.abs(phi_vals_h_vpa)
        bpar_array_vpa = np.abs(bpar_vals_h_vpa)
        apar_array_h_vpa = np.abs(100*(apar_vals_h_vpa - analytic_apar_h)/analytic_apar_h)
        apar_array_gbar_vpa = np.abs(100*(apar_vals_gbar_vpa - analytic_apar_gbar)/analytic_apar_gbar)
        apar_array_h_vperp = np.abs(100*(apar_vals_h_vperp - analytic_apar_h)/analytic_apar_h)
        apar_array_gbar_vperp = np.abs(100*(apar_vals_gbar_vperp - analytic_apar_gbar)/analytic_apar_gbar)
        ## Plot vpa-res test
        for vpamax_idx, unique_vpamax_val in enumerate(unique_vpamax_h):
            idxs = np.argwhere(np.abs(vpamax_array_h-unique_vpamax_val)<1E-4)
            ax1.scatter(nvpa_array_h[idxs], apar_array_h_vpa[idxs], marker=marker_list[vpamax_idx], s=marker_size)
        for vpamax_idx, unique_vpamax_val in enumerate(unique_vpamax_gbar):
            idxs = np.argwhere(np.abs(vpamax_array_gbar-unique_vpamax_val)<1E-4)
            ax2.scatter(nvpa_array_gbar[idxs], apar_array_gbar_vpa[idxs], marker=marker_list[vpamax_idx],
                        s=marker_size, label=r"$\tilde{v}_{\parallel, \rm max}$="+str(unique_vpamax_val))
        ## Plot vperp-res test
        for vperpmax_idx, unique_vperpmax_val in enumerate(unique_vperpmax_h):
            idxs = np.argwhere(np.abs(vperpmax_array_h-unique_vperpmax_val)<1E-4)
            ax3.scatter(nmu_array_h[idxs], apar_array_h_vperp[idxs], marker=marker_list[vperpmax_idx], s=marker_size)
        for vperpmax_idx, unique_vperpmax_val in enumerate(unique_vperpmax_gbar):
            idxs = np.argwhere(np.abs(vperpmax_array_gbar-unique_vperpmax_val)<1E-4)
            # print("unique_vperpmax_val = ", unique_vperpmax_val)
            # # print("vperpmax_idx = ", vperpmax_idx)
            # print("nmu_array_gbar[idxs] = ", nmu_array_gbar[idxs])
            # print("apar_array_gbar_vperp[idxs] = ", apar_array_gbar_vperp[idxs])
            ax4.scatter(nmu_array_gbar[idxs], apar_array_gbar_vperp[idxs], marker=marker_list[vperpmax_idx],
                        s=marker_size, label=r"$\tilde{v}_{\perp, \rm max}$="+str(unique_vperpmax_val))

        ## And now kperp
        # Need to sort the apar vals since apar depends on kperp

        apar_array_h_kperp = np.abs(100*(apar_vals_h_kperp - analytic_apar_h_kperp)/analytic_apar_h_kperp)
        apar_array_gbar_kperp = np.abs(100*(apar_vals_gbar_kperp - analytic_apar_gbar_kperp)/analytic_apar_gbar_kperp)
        ax5.scatter(kperp_vals_h, apar_array_h_kperp, marker="s", s=marker_size)
        ax6.scatter(kperp_vals_gbar, apar_array_gbar_kperp, marker="s", s=marker_size)

        for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.tick_params("x", which="major", length=my_xticklength, width=my_xtickwidth, direction="in")
            ax.tick_params("x", which="minor", length=my_xminorticklength, width=my_xminortickwidth, direction="in")
            ax.tick_params("y", which="major", length=my_yticklength, width=my_ytickwidth, direction="in")
            ax.tick_params("y", which="minor", length=my_yminorticklength, width=my_yminortickwidth, direction="in")

        for ax in [ax1, ax2, ax3, ax4]:
            ax.set_xticks([10, 100])
        for ax in [ax1, ax3]:
            ax.set_xticklabels([])
        for ax in [ax2, ax4]:
            ax.set_xticklabels([r"$10^1$", r"$10^2$"], fontsize=xticklabel_fontsize)

        for ax in [ax5, ax6]:
            ax.set_xticks([1e-4, 1e-2, 1, 100])
            ax.set_xticks([1e-3, 1e-1, 10], minor=True)
            ax.set_xticklabels([], minor=True)

        ax5.set_xticklabels([])
        ax6.set_xticklabels([r"$10^{-4}$", r"$10^{-2}$", r"$1$", r"$10^{2}$"], fontsize=xticklabel_fontsize)

        ax2.set_ylim((1E-18, 3))
        ax4.set_ylim((1E-19, 2))

        ax1.set_yticks([1e-8, 1e-4, 1])
        ax1.set_yticklabels([r"$10^{-8}$",
            r"$10^{-4}$", r"$1$"], fontsize=yticklabel_fontsize)
        ax1.set_yticks([1e-6, 1e-2, 100], minor=True)
        ax1.set_yticklabels([], minor=True)

        ax2.set_yticks([1e-16, 1e-12, 1e-8, 1e-4, 1])
        ax2.set_yticklabels([r"$10^{-16}$", r"$10^{-12}$", r"$10^{-8}$",
            r"$10^{-4}$", r"$1$"], fontsize=yticklabel_fontsize)
        ax2.set_yticks([1e-14, 1e-10, 1e-6, 1e-2], minor=True)
        ax2.set_yticklabels([], minor=True)

        ax3.set_yticks([1e-8, 1e-4, 1])
        ax3.set_yticklabels([r"$10^{-8}$",
            r"$10^{-4}$", r"$1$"], fontsize=yticklabel_fontsize)
        ax3.set_yticks([1e-6, 1e-2, 100], minor=True)
        ax3.set_yticklabels([], minor=True)

        ax4.set_yticks([1e-16, 1e-12, 1e-8, 1e-4, 1])
        ax4.set_yticklabels([r"$10^{-16}$", r"$10^{-12}$", r"$10^{-8}$",
            r"$10^{-4}$", r"$1$"], fontsize=yticklabel_fontsize)
        ax4.set_yticks([1e-18, 1e-14, 1e-10, 1e-6, 1e-2], minor=True)
        ax4.set_yticklabels([], minor=True)

        ax5.set_yticks([1e-8, 1e-6, 1e-4])
        ax5.set_yticklabels([r"$10^{-8}$", r"$10^{-6}$",
            r"$10^{-4}$"], fontsize=yticklabel_fontsize)
        # ax5.set_yticks([1e-9, 1e-7, 1e-5, 1e-3], minor=True)
        # ax5.set_yticklabels([], minor=True)

        ax6.set_yticks([1e-12, 1e-11, 1e-10, 1e-9])
        ax6.set_yticklabels([r"$10^{-12}$", r"$10^{-11}$", r"$10^{-10}$",
            r"$10^{-9}$"], fontsize=yticklabel_fontsize)
        # ax6.set_yticks([1e-18, 1e-14, 1e-10, 1e-6, 1e-2], minor=True)
        # ax6.set_yticklabels([], minor=True)

        for ax in [ax2, ax4]:
            ax.legend(loc="best", ncol=2)

        for ax in [ax2]:
            ax.set_xlabel(r"$n_{v \parallel}$", fontsize=xlabel_fontsize)
        for ax in [ax4]:
            ax.set_xlabel(r"$n_{v \perp}$", fontsize=xlabel_fontsize)
        for ax in [ax6]:
            ax.set_xlabel(r"$\tilde{k}_\perp$", fontsize=xlabel_fontsize)

        ax1.set_ylabel(r"$\Delta \tilde{A}_{\parallel k} (\%)$ ", fontsize=ylabel_fontsize)
        ax2.set_ylabel(r"$\Delta \tilde{A}_{\parallel k} (\%)$ ", fontsize=ylabel_fontsize)
        plt.savefig("apar_field_solve.eps")
        plt.close()

        return

    analytic_phi_h, analytic_apar_h, analytic_bpar_h = calculate_fields_zvpajzeroexp(1, 1, "h")
    analytic_phi_gbar, analytic_apar_gbar, analytic_bpar_gbar = calculate_fields_zvpajzeroexp(1, 1, "gbar")
    kperp_vals = np.logspace(-4, 2, 20)

    analytic_apar_h_kperp = []
    analytic_apar_gbar_kperp = []
    for kperp_val in kperp_vals:
        phi_h, apar_h, bpar_h = calculate_fields_zvpajzeroexp(kperp_val, 1, "h")
        phi_gbar, apar_gbar, bpar_gbar = calculate_fields_zvpajzeroexp(kperp_val, 1, "gbar")
        analytic_apar_h_kperp.append(apar_h)
        analytic_apar_gbar_kperp.append(apar_gbar)
    analytic_apar_h_kperp = np.array(analytic_apar_h_kperp)
    analytic_apar_gbar_kperp = np.array(analytic_apar_gbar_kperp)
    make_analytic_apar_plot()
    # print("solution with h. apar = ", analytic_apar_h)
    # print("solution with gbar. apar = ", analytic_apar_gbar)
    # print("solution with h. analytic_apar_h_kperp = ", analytic_apar_h_kperp)
    # print("solution with gbar. analytic_apar_gbar_kperp = ", analytic_apar_gbar_kperp)


    vpamax_vals_h, nvpa_vals_h, phi_vals_h_vpa, apar_vals_h_vpa, bpar_vals_h_vpa = get_fields_from_outnc_files(make_sims.apar_h_vpa_scan_folder)
    vpamax_vals_gbar, nvpa_vals_gbar, phi_vals_gbar_vpa, apar_vals_gbar_vpa, bpar_vals_gbar_vpa = get_fields_from_outnc_files(make_sims.apar_gbar_vpa_scan_folder)
    vperpmax_vals_h, nmu_vals_h, phi_vals_h_vperp, apar_vals_h_vperp, bpar_vals_h_vperp = get_fields_from_outnc_files(make_sims.apar_h_vperp_scan_folder)
    vperpmax_vals_gbar, nmu_vals_gbar, phi_vals_gbar_vperp, apar_vals_gbar_vperp, bpar_vals_gbar_vperp = get_fields_from_outnc_files(make_sims.apar_gbar_vperp_scan_folder)
    kperp_vals_h, dum, phi_vals_h_kperp, apar_vals_h_kperp, bpar_vals_h_kperp = get_fields_from_outnc_files(make_sims.apar_h_kperp_scan_folder, kperp=True)
    kperp_vals_gbar, dum, phi_vals_gbar_kperp, apar_vals_gbar_kperp, bpar_vals_gbar_kperp = get_fields_from_outnc_files(make_sims.apar_gbar_kperp_scan_folder, kperp=True)
    print("kperp_vals_gbar = ", kperp_vals_gbar)
    print("max(abs(phi_vals_h_kperp)) = ", np.max(abs(phi_vals_h_kperp)))
    print("max(abs(bpar_vals_h_kperp)) = ", np.max(abs(bpar_vals_h_kperp)))
    print("max(abs(phi_vals_gbar_kperp)) = ", np.max(abs(phi_vals_gbar_kperp)))
    arg = np.argmax(abs(phi_vals_gbar_kperp))
    print("kperp of max phi = ", kperp_vals_gbar[arg])
    print("max(abs(bpar_vals_gbar_kperp)) = ", np.max(abs(bpar_vals_gbar_kperp)))
    sort_idxs = np.argsort(kperp_vals_h)
    kperp_vals_h = kperp_vals_h[sort_idxs]
    apar_vals_h_kperp = apar_vals_h_kperp[sort_idxs]
    sort_idxs = np.argsort(kperp_vals_gbar)
    kperp_vals_gbar = kperp_vals_gbar[sort_idxs]
    apar_vals_gbar_kperp = apar_vals_gbar_kperp[sort_idxs]

    # print("apar_vals_h_kperp = ", apar_vals_h_kperp)
    # print("apar_vals_gbar_kperp = ", apar_vals_gbar_kperp)
    # print("phi_vals_h = ", phi_vals_h)
    # print("bpar_vals_h = ", bpar_vals_h)
    # print("phi_vals_gbar = ", phi_vals_gbar)
    # print("apar_vals_gbar = ", apar_vals_gbar)
    # print("bpar_vals_gbar = ", bpar_vals_gbar)
    make_plot()

    # make_plot("vperp")
    return

if __name__ == "__main__":
    print("Hello world")
    vspace_res_test_field_solve_for_thesis()
    test_field_solve_apar_for_thesis()
