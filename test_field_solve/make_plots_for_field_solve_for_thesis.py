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

    y_ticklabelfontsize = 25
    x_ticklabelfontsize = y_ticklabelfontsize
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
        ax1.set_yticklabels(ax1_yticklabels, fontsize=y_ticklabelfontsize)
    if ax2_yticks is not None:
        ax2.set_yticks(ax2_yticks)
        ax2.set_yticklabels(ax2_yticklabels, fontsize=y_ticklabelfontsize)
    ax1.tick_params("x", length=my_xticklength, width=my_xtickwidth, direction="out")
    ax2.tick_params("x", length=my_xticklength, width=my_xtickwidth, direction="out")
    ax1.tick_params("y", length=my_yticklength, width=my_ytickwidth, direction="out")
    ax1.set_xticks([-1, 0, 1])
    ax1.set_xticklabels([])
    ax2.set_xticks([-1, 0, 1])
    ax2.set_xticklabels(["-1", "0", "1"], fontsize=x_ticklabelfontsize)
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

def make_field_solve_test_plots_for_thesis():
    """Make plots for thesis:
    fapar=fbpar=0: phi benchmark
    fapar=1, fbpar=0: benchmarking phi, apar (where apar!=0)
    fapar=0, fbpar=1: benchmarking phi, bpar
    fapar=fbpar=1: benchmarking phi, apar, bpar (where a)apar=0, b)apar!=0 ) """

    compare_phi_for_thesis([gs2_sim_fapar0_fbpar0, gs2_sim_fapar0_fbpar0_lr, gs2_sim_fapar0_fbpar0_vlr,
                 stella_sim_fapar0_fbpar0, stella_sim_fapar0_fbpar0_lr, stella_sim_fapar0_fbpar0_vlr],
                ["gs2", "gs2", "gs2",
                 "stella", "stella", "stella"],
                ["gs2 (nguass=12, negrid=16)", "gs2 (nguass=6, negrid=8)",
                 "gs2 (nguass=3, negrid=4)",
                 "stella (nvgrid=24, nmu=12)",
                 "stella (nvgrid=12, nmu=6)",
                 "stella (nvgrid=6, nmu=3)",
                 ],
                 save_name="field_solve_test_fapar0_fbpar0.eps",
                 ax1_yticks=[-4, 0, 4],
                 ax1_yticklabels=[r"$-4$", r"$0$", r"$4$"],
                 ax2_yticks=[10**(-4), 10**(-3), 10**(-2), 10**(-1), 10**(0), 10**(1)],
                 ax2_yticklabels=[r"$10^{-4}$", r"$10^{-3}$", r"$10^{-2}$", r"$10^{-1}$",
                                     r"$10^{0}$", r"$10^{1}$"],
         )
    # ### fapar=1, fbpar=0
    # compare_apar([gs2_sim_fapar1_fbpar0, stella_sim_fapar1_fbpar0],
    #             ["gs2", "stella"],
    #             ["gs2", "stella"], title="fapar=1, fbpar=0")
    # ### fapar=0, fbpar=1
    compare_phi_for_thesis([gs2_sim_fapar0_fbpar1, stella_sim_fapar0_fbpar1],
                ["gs2", "stella"],
                ["gs2", "stella"],
                save_name="field_solve_test_fapar0_fbpar1_phi.eps",)
    compare_bpar_for_thesis([gs2_sim_fapar0_fbpar1, stella_sim_fapar0_fbpar1],
                ["gs2", "stella"],
                ["gs2", "stella"],
                save_name="field_solve_test_fapar0_fbpar1_bpar.eps",)
    ### fapar=1, fbpar=1
    compare_phi_for_thesis([gs2_sim_fapar1_fbpar1, gs2_sim_fapar1_fbpar1_lr, gs2_sim_fapar1_fbpar1_vlr,
                stella_sim_fapar1_fbpar1, stella_sim_fapar1_fbpar1_lr, stella_sim_fapar1_fbpar1_vlr
                ],
                ["gs2", "gs2", "gs2",
                 "stella", "stella", "stella"],
                ["gs2 (nguass=12, negrid=16)", "gs2 (nguass=6, negrid=8)",
                 "gs2 (nguass=3, negrid=4)",
                 "stella (nvgrid=24, nmu=12)",
                 "stella (nvgrid=12, nmu=6)",
                 "stella (nvgrid=6, nmu=3)",
                 ]
                 , save_name="field_solve_test_fapar1_fbpar1_phi.eps",)
    # compare_apar([gs2_sim_fapar1_fbpar1, gs2_sim_fapar1_fbpar1_lr, gs2_sim_fapar1_fbpar1_vlr,
    #             stella_sim_fapar1_fbpar1, stella_sim_fapar1_fbpar1_lr, stella_sim_fapar1_fbpar1_vlr
    #             ],
    #             ["gs2", "gs2", "gs2",
    #              "stella", "stella", "stella"],
    #             ["gs2 (nguass=12, negrid=16)", "gs2 (nguass=6, negrid=8)",
    #              "gs2 (nguass=3, negrid=4)",
    #              "stella (nvgrid=24, nmu=12)",
    #              "stella (nvgrid=12, nmu=6)",
    #              "stella (nvgrid=6, nmu=3)",
    #              ]
    #              , title="fapar=1, fbpar=1")
    compare_bpar_for_thesis([gs2_sim_fapar1_fbpar1, gs2_sim_fapar1_fbpar1_lr, gs2_sim_fapar1_fbpar1_vlr,
                stella_sim_fapar1_fbpar1, stella_sim_fapar1_fbpar1_lr, stella_sim_fapar1_fbpar1_vlr
                ],
                ["gs2", "gs2", "gs2",
                 "stella", "stella", "stella"],
                ["gs2 (nguass=12, negrid=16)", "gs2 (nguass=6, negrid=8)",
                 "gs2 (nguass=3, negrid=4)",
                 "stella (nvgrid=24, nmu=12)",
                 "stella (nvgrid=12, nmu=6)",
                 "stella (nvgrid=6, nmu=3)",
                 ]
                 , save_name="field_solve_test_fapar1_fbpar1_bpar.eps",)
    return

def test_field_solve_in_h():
    """ """
    compare_phi_for_thesis([
                stella_sim_fapar1_fbpar1, stella_sim_fapar1_fbpar1_h
                            ],
                ["stella", "stella"],
                [
                 "stella (nvgrid=24, nmu=12)",
                 "stella (nvgrid=12, nmu=6) (in h)",
                 ]
                 , save_name="h_field_solve_test_fapar1_fbpar1_phi.eps",)
    compare_bpar_for_thesis([
                stella_sim_fapar1_fbpar1, stella_sim_fapar1_fbpar1_h
                            ],
                ["stella", "stella"],
                [
                 "stella (nvgrid=24, nmu=12)",
                 "stella (nvgrid=12, nmu=6) (in h)",
                 ]
                 , save_name="h_field_solve_test_fapar1_fbpar1_bpar.eps",)
    return

def test_analytic_field_solve_in_h():
    """ """
    m_ion = 1
    m_electron = 2.8E-4
    B = 1

    def calculate_phi(kperp):
        """ """
        b_ion = 0.5*kperp*kperp*m_ion/(B*B)
        b_electron = 0.5*kperp*kperp*m_electron/(B*B)
        gamzero_ion = np.exp(-b_ion)*iv(0,b_ion)
        gamzero_electron = np.exp(-b_electron)*iv(0,b_electron)

        phi = 0.5*(gamzero_ion + gamzero_electron)

        return phi

    def calculate_bpar(kperp, beta):
        """ """
        b_ion = 0.5*kperp*kperp*m_ion/(B*B)
        b_electron = 0.5*kperp*kperp*m_electron/(B*B)
        gamone_ion = np.exp(-b_ion)*(iv(0,b_ion) - iv(1,b_ion) )
        gamone_electron = np.exp(-b_electron)*(iv(0,b_electron) - iv(1,b_electron) )
        bpar = (-0.5*beta/B) * (gamone_ion - gamone_electron )
        return bpar

    # outnc_longname = "sims/stella_slab_fapar1_fbpar1_h/input_nvpa48_nmu24.out.nc"
    outnc_longname = "sims/stella_slab_fapar1_fbpar1_h/input_vpamax4_vperpmax4.out.nc"
    field_key_stella = "phi_vs_t"
    (z, field_vs_t) = extract_data_from_ncdf_with_xarray(outnc_longname,
            'zed', field_key_stella)
    field_t0 = np.array(field_vs_t[0,0,:,0,0,:])
    # shape is now (zed, ri)
    print("stella phi = ", field_t0)
    phi = calculate_phi(1)
    print("analytic phi = ", phi)
    print("stella/analytic = ", field_t0[0]/phi)

    field_key_stella = "bpar_vs_t"
    (z, field_vs_t) = extract_data_from_ncdf_with_xarray(outnc_longname,
            'zed', field_key_stella)
    field_t0 = np.array(field_vs_t[0,0,:,0,0,:])
    print("stella bpar = ", field_t0)
    kperp = 1; beta=1
    bpar = calculate_bpar(kperp, beta)
    print("analytic bpar = ", bpar)
    print("stella/analytic = ", field_t0[0]/bpar)

    return

def vspace_res_test_field_solve_in_h_for_thesis():
    """Plot % error in phi, bpar and abs(apar) in various cases:
    (1) Fixed vperpmax, nmu. plotting (error in field) vs nvpa for several vpamax
    (2) Fixed vpamax, nvpa, plotting (error in field) vs nvmu for several vperpmax
    For each of these, plot the result for field solve in h and field solve in h.
    """

    def make_plot(which):
        """ """
        marker_size = 12
        label_fontsize = 40
        legend_fontsize = 14
        marker_list = ["s", "o", "P", "X", "v", "^", "<", ">", "1", "2", "3"]
        lw_list = [4, 3, 3, 2]
        ls_list = ["-", "--", "-.", (0, (4,1,2,1))]

        my_xticklength = 7
        my_xtickwidth = 2
        my_xticklength_minor = 4
        my_xtickwidth_minor = 1

        x_ticklabelfontsize = 20
        y_ticklabelfontsize = 20
        y_labelfontsize = 30
        x_labelfontsize = 30
        bracket_fontsize = 70
        itg_fontsize = 30

        top = 0.98
        left = 0.14
        right = 0.98
        bottom = 0.13
        vspace = 0.02
        hspace = 0.04
        height = (top - bottom - 2*vspace)/3
        row3_bottom = bottom
        row2_bottom = bottom + height + vspace
        row1_bottom = row2_bottom + height + vspace
        width = (right - left - hspace)/2
        col1_left = left
        col2_left = left + width + hspace

        fig = plt.figure(figsize=(12,16))
        ax1 = fig.add_axes((col1_left, row1_bottom, width, height))
        ax2 = fig.add_axes((col2_left, row1_bottom, width, height))
        ax3 = fig.add_axes((col1_left, row2_bottom, width, height))
        ax4 = fig.add_axes((col2_left, row2_bottom, width, height))
        ax5 = fig.add_axes((col1_left, row3_bottom, width, height))
        ax6 = fig.add_axes((col2_left, row3_bottom, width, height))

        if which=="vpa":
            for ax in [ax5, ax6]:
                ax.set_xlabel(r"nvpa", fontsize=x_labelfontsize)
        elif which=="vperp":
            for ax in [ax5, ax6]:
                ax.set_xlabel(r"nvperp", fontsize=x_labelfontsize)

        ax1.set_ylabel(r"$\Delta \tilde{\phi}_k (\%)$ ", fontsize=y_labelfontsize)
        ax3.set_ylabel(r"$\Delta \tilde{A}_{\parallel k}$ ", fontsize=y_labelfontsize)
        ax5.set_ylabel(r"$\Delta \tilde{B}_{\parallel k} (\%)$ ", fontsize=y_labelfontsize)
        if which=="vpa":
            plt.savefig("vpa_res_field_solve.eps")
        if which=="vperp":
            plt.savefig("vperp_res_field_solve.eps")
        plt.close()

    make_plot("vpa")
    make_plot("vperp")
    return


if __name__ == "__main__":
    print("Hello world")
    # make_field_solve_test_plots_for_thesis()
    # test_field_solve_in_h()
    # test_analytic_field_solve_in_h()
    vspace_res_test_field_solve_in_h_for_thesis()
