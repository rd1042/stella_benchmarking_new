import sys
import pickle
import numpy as np
import glob
import re
import matplotlib.pyplot as plt
from itertools import cycle
import xarray as xr
from scipy.interpolate import interp1d

sys.path.append("../postprocessing_tools")
import helper_linear_sims as help_lin
from extract_sim_data import get_omega_data, get_phiz_data, get_aparz_data, get_bparz_data
from plotting_helper import plot_phi_z_for_sim, plot_apar_z_for_sim, plot_bpar_z_for_sim
from helper_ncdf_new import view_ncdf_variables, extract_data_from_ncdf_with_xarray
import plot_2d_utils as plot2dutils

beta_04_for_thesis_folder = "sims/beta_0.04000_for_thesis/"
beta_0_for_thesis_folder = "sims/beta_0_for_thesis/"
default_cmap = plt.get_cmap("tab10")

def make_comparison_beta004(diff=False):
    """ """
    marker_size = 12
    legend_fontsize = 24
    marker_list = ["s", "o", "P", "X", "v", "^", "<", ">", "1", "2", "3"]
    lw_list = [6, 4, 2]
    ls_list = ["-", "--", "-.", (0, (4,1,2,1))]

    my_xticklength = 7
    my_xtickwidth = 2
    my_xticklength_minor = 4
    my_xtickwidth_minor = 1
    my_yticklength = 7
    my_ytickwidth = 2
    my_yticklength_minor = 4
    my_ytickwidth_minor = 1

    x_ticklabelfontsize = 20
    y_ticklabelfontsize = 20
    y_labelfontsize = 30
    x_labelfontsize = 30
    bracket_fontsize = 70
    itg_fontsize = 30

    top = 0.98
    left = 0.14
    right = 0.98
    bottom = 0.1
    vspace = 0.03
    height = (top - bottom - 2*vspace)/3
    row2_bottom = bottom + height + vspace
    row1_bottom = row2_bottom + height + vspace
    width = right - left
    outnc_longnames= [
        beta_04_for_thesis_folder + "gs2_nt64_np4.out.nc",
        beta_04_for_thesis_folder + "stella_implicit_np4_nt64.out.nc",
        beta_04_for_thesis_folder + "stella_explicit_zupw0_src_h_fapar1_fbpar1.out.nc",
        #beta_04_for_thesis_folder + "stella_implicit_np4_nt64_tupw0.out.nc",
        ]
    labels =[
             r"GS2 ($n_\theta=64$, $n_{period}=4$)",
              r"stella (implicit) ($n_\theta=64$, $n_{period}=4$)",
              r"stella (explicit) ($n_\theta=64$, $n_{period}=4$)",
              #r"stella (implicit) ($n_\theta=64$, $n_{period}=4$, $u_t=0$)",
                ]
    col_list = [default_cmap(1), default_cmap(3), default_cmap(4),
                # default_cmap(4)
                ]
    sim_types = ["gs2", "stella", "stella",
                 # "stella"
                 ]

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
    bmag_list = []
    bmag_z_list = []
    for sim_idx, outnc_longname in enumerate(outnc_longnames):
        sim_type = sim_types[sim_idx]
        ## For each sim, we want to compare:
        # Mode structures
        # Omega(t) (are we converged?)
        # Omega (print to screen?)
        sim_shortname = re.split("/", outnc_longname)[-1]
        sim_shortname = re.split(".out.nc", sim_shortname)[0]
        sim_longname = re.split(".out.nc", outnc_longname)[0]
        # view_ncdf_variables(outnc_longname)

        z, real_phi, imag_phi = get_phiz_data(sim_longname, sim_type)
        abs_phi = help_lin.get_abs(real_phi, imag_phi)
        z, real_apar, imag_apar = get_aparz_data(sim_longname, sim_type)
        abs_apar = help_lin.get_abs(real_apar, imag_apar)
        z, real_bpar, imag_bpar = get_bparz_data(sim_longname, sim_type)
        abs_bpar = help_lin.get_abs(real_bpar, imag_bpar)
        max_abs_phi = np.max(abs_phi)

        sim_shortnames.append(sim_shortname)
        zvals_list.append(z)
        phiz_list.append(abs_phi/max_abs_phi)
        aparz_list.append(abs_apar/max_abs_phi)
        bparz_list.append(abs_bpar/max_abs_phi)

    z_ref = zvals_list[0]
    phiz_ref = phiz_list[0]
    aparz_ref = aparz_list[0]
    bparz_ref = bparz_list[0]

    fig = plt.figure(figsize=(12,12))
    ax1 = fig.add_axes((left, row1_bottom, width, height))
    ax2 = fig.add_axes((left, row2_bottom, width, height))
    ax3 = fig.add_axes((left, bottom, width, height))

    for sim_idx, sim_shortname in enumerate(sim_shortnames):
        if not diff:
            ax1.plot(zvals_list[sim_idx]/np.pi, phiz_list[sim_idx], label=labels[sim_idx],
                c=col_list[sim_idx], ls=ls_list[sim_idx],lw=lw_list[sim_idx])
            ax2.plot(zvals_list[sim_idx]/np.pi, aparz_list[sim_idx], label=labels[sim_idx],
                c=col_list[sim_idx], ls=ls_list[sim_idx],lw=lw_list[sim_idx])
            ax3.plot(zvals_list[sim_idx]/np.pi, bparz_list[sim_idx], label=labels[sim_idx],
                c=col_list[sim_idx], ls=ls_list[sim_idx],lw=lw_list[sim_idx])
        else:
            if sim_idx!=0:
                ## Whether or not this is the reference beta scan, plot Omega-Omega_ref
                # Find where beta vals match.
                # print("phiz_list[sim_idx] = ", phiz_list[sim_idx] )
                # print("phiz_ref = ", phiz_ref )
                z = zvals_list[sim_idx]
                f_phi = interp1d(np.array(z_ref), np.array(phiz_ref))
                phi_interp = f_phi(np.array(z))
                f_apar = interp1d(z_ref, aparz_ref)
                apar_interp = f_apar(np.array(z))
                # z0_idx = int(len(z_ref)*0.5)
                # print("z_ref[z0_idx] = ", z_ref[z0_idx])
                # print("z[z0_idx] = ", z[z0_idx])
                # apar_interp[z0_idx] = aparz_ref[z0_idx]
                f_bpar = interp1d(z_ref, bparz_ref)
                bpar_interp = f_bpar(np.array(z))
                # print("z = ", z)
                # print("phiz_list[sim_idx] = ", phiz_list[sim_idx])
                phiz_diff = (phiz_list[sim_idx] - phi_interp) / phi_interp
                aparz_diff = (aparz_list[sim_idx] - apar_interp) / apar_interp
                bparz_diff = (bparz_list[sim_idx] - bpar_interp) / bpar_interp
                ax1.plot(z/np.pi, phiz_diff*100, label=labels[sim_idx],
                         lw=lw_list[sim_idx], ls=ls_list[sim_idx], markersize=marker_size,
                         c=col_list[sim_idx])
                ax2.plot(z/np.pi, aparz_diff*100, label=labels[sim_idx],
                         lw=lw_list[sim_idx], ls=ls_list[sim_idx], markersize=marker_size,
                         c=col_list[sim_idx])
                ax3.plot(z/np.pi, bparz_diff*100, label=labels[sim_idx],
                         lw=lw_list[sim_idx], ls=ls_list[sim_idx], markersize=marker_size,
                         c=col_list[sim_idx])

    for ax in [ax1, ax2, ax3]:
        ax.set_xlim([-7, 7])
        ax.tick_params("x", which="major", length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.tick_params("x", which="minor", length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in")
        ax.tick_params("y", which="major", length=my_yticklength, width=my_ytickwidth, direction="in")
        ax.tick_params("y", which="minor", length=my_yticklength_minor, width=my_ytickwidth_minor, direction="in")
        ax.set_xticks([-6, -4, -2, 0, 2, 4, 6])
        ax.set_xticks([-7, -5, -3, -1, 1, 3, 5, 7], minor=True)
        ax.set_xticklabels([], minor=True)
    ax2.set_xticklabels([], fontsize=x_ticklabelfontsize)
    ax1.set_xticklabels([], fontsize=x_ticklabelfontsize)
    ax3.set_xticklabels([r"$-6\pi$", r"$-4\pi$", r"$-2\pi$", r"$0$", r"$2\pi$", r"$4\pi$", r"$6\pi$"], fontsize=x_ticklabelfontsize)

    ax1.legend(loc="best", fontsize=legend_fontsize)
    if diff:
        ax1.set_ylim([-40, 58])
        ax2.set_ylim([-20, 5])
        ax3.set_ylim([-35, 32])
        ax1.set_yticks([-40, 0, 40,])
        ax1.set_yticklabels([r"$-40$", r"$0$", r"$40$"], fontsize=y_ticklabelfontsize)
        ax1.set_yticks([-20, 20], minor=True)
        ax1.set_yticklabels([], minor=True)
        ax2.set_yticks([-20, -10, 0, 10,])
        ax2.set_yticklabels([r"$-20$", r"$-10$", r"$0$", r"$10$"], fontsize=y_ticklabelfontsize)
        ax2.set_yticks([-15, -5, 5], minor=True)
        ax2.set_yticklabels([], minor=True)
        ax3.set_yticks([-20, 0, 20,])
        ax3.set_yticklabels([r"$-20$", r"$0$", r"$20$"], fontsize=y_ticklabelfontsize)
        ax3.set_yticks([-30, -10, 10, 30], minor=True)
        ax1.set_ylabel(r"$\Delta \vert \tilde{\phi}_k \vert$ ($\%$)", fontsize=y_labelfontsize)
        ax2.set_ylabel(r"$\Delta \vert \tilde{A}_{\parallel k} \vert$ ($\%$)", fontsize=y_labelfontsize)
        ax3.set_ylabel(r"$\Delta \vert \tilde{B}_{\parallel k} \vert$ ($\%$)", fontsize=y_labelfontsize)
    else:
        ax1.set_ylim([-0.02, 1.02])
        ax2.set_ylim([-0.003, 0.068])
        ax3.set_ylim([-0.001, 0.026])
        ax1.set_ylabel(r"$\vert \tilde{\phi}_k \vert$", fontsize=y_labelfontsize)
        ax2.set_ylabel(r"$\vert \tilde{A}_{\parallel k} \vert$", fontsize=y_labelfontsize)
        ax3.set_ylabel(r"$\vert \tilde{B}_{\parallel k} \vert$", fontsize=y_labelfontsize)
        ax1.set_yticks([0, 0.5, 1,])
        ax1.set_yticklabels([r"$0$", r"$0.5$", r"$1$"], fontsize=y_ticklabelfontsize)
        ax1.set_yticks([0.25, 0.75], minor=True)
        ax1.set_yticklabels([], minor=True)
        ax2.set_yticks([0, 0.03, 0.06,])
        ax2.set_yticklabels([r"$0$", r"$0.03$", r"$0.06$"], fontsize=y_ticklabelfontsize)
        ax2.set_yticks([0.01, 0.02, 0.04, 0.05], minor=True)
        ax2.set_yticklabels([], minor=True)
        ax3.set_yticks([0, 0.01, 0.02,])
        ax3.set_yticklabels([r"$0$", r"$0.01$", r"$0.02$"], fontsize=y_ticklabelfontsize)
        ax3.set_yticks([0.005, 0.015], minor=True)
        ax3.set_yticklabels([], minor=True)

    ax3.set_xlabel(r"$z$", fontsize=x_labelfontsize)


    for ax in [ax1, ax2, ax3]:
        ax.grid(True)

    if diff:
        plt.savefig("beta_004_mode_comparison_diff.eps")
    else:
        plt.savefig("beta_004_mode_comparison_fields.eps")
    return

def make_comparison_beta0(diff=False):
    """ """
    marker_size = 12
    legend_fontsize = 16
    marker_list = ["s", "o", "P", "X", "v", "^", "<", ">", "1", "2", "3"]
    lw_list = [6, 4, 2]
    ls_list = ["-", "--", "-.", (0, (4,1,2,1))]

    my_xticklength = 7
    my_xtickwidth = 2
    my_xticklength_minor = 4
    my_xtickwidth_minor = 1
    my_yticklength = 7
    my_ytickwidth = 2
    my_yticklength_minor = 4
    my_ytickwidth_minor = 1

    x_ticklabelfontsize = 20
    y_ticklabelfontsize = 20
    y_labelfontsize = 30
    x_labelfontsize = 30
    bracket_fontsize = 70
    itg_fontsize = 30

    top = 0.98
    left = 0.14
    right = 0.98
    bottom = 0.1
    vspace = 0.03
    height = (top - bottom)
    width = right - left

    outnc_longnames= [
        beta_0_for_thesis_folder + "gs2_nt64_np4.out.nc",
        beta_0_for_thesis_folder + "stella_implicit_np4_nt64.out.nc",
        beta_0_for_thesis_folder + "stella_explicit_zupw0_src_h.out.nc",
        ]
    labels =[
             r"GS2 ($n_\theta=64$, $n_{period}=4$)",
              r"stella (implicit) ($n_\theta=64$, $n_{period}=4$)",
              r"stella (explicit) ($n_\theta=64$, $n_{period}=4$)",
                ]
    sim_types = ["gs2",
                 "stella",
                 "stella"
                 ]
    col_list = [default_cmap(1), default_cmap(3), default_cmap(4)]

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
    bmag_list = []
    bmag_z_list = []
    for sim_idx, outnc_longname in enumerate(outnc_longnames):
        sim_type = sim_types[sim_idx]
        ## For each sim, we want to compare:
        # Mode structures
        # Omega(t) (are we converged?)
        # Omega (print to screen?)
        sim_shortname = re.split("/", outnc_longname)[-1]
        sim_shortname = re.split(".out.nc", sim_shortname)[0]
        sim_longname = re.split(".out.nc", outnc_longname)[0]
        # view_ncdf_variables(outnc_longname)

        z, real_phi, imag_phi = get_phiz_data(sim_longname, sim_type)
        abs_phi = help_lin.get_abs(real_phi, imag_phi)
        max_abs_phi = np.max(abs_phi)

        sim_shortnames.append(sim_shortname)
        zvals_list.append(z)
        phiz_list.append(abs_phi/max_abs_phi)

    z_ref = zvals_list[0]
    phiz_ref = phiz_list[0]

    fig = plt.figure(figsize=(12,8))
    ax1 = fig.add_axes((left, bottom, width, height))

    for sim_idx, sim_shortname in enumerate(sim_shortnames):
        if not diff:
            ax1.plot(zvals_list[sim_idx]/np.pi, phiz_list[sim_idx], label=labels[sim_idx],
                c=col_list[sim_idx], ls=ls_list[sim_idx],lw=lw_list[sim_idx])
        else:
            if sim_idx!=0:
                ## Whether or not this is the reference beta scan, plot Omega-Omega_ref
                # Find where beta vals match.
                # print("phiz_list[sim_idx] = ", phiz_list[sim_idx] )
                # print("phiz_ref = ", phiz_ref )
                z = zvals_list[sim_idx]
                f_phi = interp1d(np.array(z_ref), np.array(phiz_ref))
                phi_interp = f_phi(np.array(z))
                phiz_diff = (phiz_list[sim_idx] - phi_interp) / phi_interp
                ax1.plot(z/np.pi, phiz_diff*100, label=labels[sim_idx],
                         lw=lw_list[sim_idx], ls=ls_list[sim_idx], markersize=marker_size,
                         c=col_list[sim_idx])

    for ax in [ax1]:
        ax.set_xlim([-7, 7])
        ax.tick_params("x", which="major", length=my_xticklength, width=my_xtickwidth, direction="in")
        ax.tick_params("x", which="minor", length=my_xticklength_minor, width=my_xtickwidth_minor, direction="in")
        ax.tick_params("y", which="major", length=my_yticklength, width=my_ytickwidth, direction="in")
        ax.tick_params("y", which="minor", length=my_yticklength_minor, width=my_ytickwidth_minor, direction="in")
        ax.set_xticks([-6, -4, -2, 0, 2, 4, 6])
        ax.set_xticks([-7, -5, -3, -1, 1, 3, 5, 7], minor=True)
        ax.set_xticklabels([], minor=True)
    ax1.set_xticklabels([r"$-6\pi$", r"$-4\pi$", r"$-2\pi$", r"$0$", r"$2\pi$", r"$4\pi$", r"$6\pi$"], fontsize=x_ticklabelfontsize)

    ax1.legend(loc="upper right", fontsize=legend_fontsize)
    if diff:
        ax1.set_ylim([-80, 79])
        ax1.set_yticks([-80, -40, 0, 40,])
        ax1.set_yticklabels([r"$-80$", r"$-40$", r"$0$", r"$40$"], fontsize=y_ticklabelfontsize)
        ax1.set_yticks([-60, -20, 20, 60], minor=True)
        ax1.set_yticklabels([], minor=True)
        ax1.set_ylabel(r"$\Delta \vert \tilde{\phi}_k \vert$ ($\%$)", fontsize=y_labelfontsize)
    else:
        ax1.set_ylim([-0.02, 1.1])
        ax1.set_ylabel(r"$\vert \tilde{\phi}_k \vert$", fontsize=y_labelfontsize)
        ax1.set_yticks([0, 0.5, 1,])
        ax1.set_yticklabels([r"$0$", r"$0.5$", r"$1$"], fontsize=y_ticklabelfontsize)
        ax1.set_yticks([0.25, 0.75], minor=True)
        ax1.set_yticklabels([], minor=True)

    ax1.set_xlabel(r"$z$", fontsize=x_labelfontsize)


    for ax in [ax1]:
        ax.grid(True)

    if diff:
        plt.savefig("beta_0_mode_comparison_diff.eps")
    else:
        plt.savefig("beta_0_mode_comparison_fields.eps")
    return

def for_thesis_beta004_benchmark_plot():
    """Compare Omega(t) and fields for simulation with nz=64, nperiod=4 for the following:
    (1) GS2
    (2) stella implicit streaming + mirror
    (3) stella explicit """
    make_comparison_beta004(diff=False)
    make_comparison_beta004(diff=True)

    return

def for_thesis_beta0_benchmark_plot():
    """Compare Omega(t) and fields for simulation with nz=64, nperiod=4 for the following:
    (1) GS2
    (2) stella implicit streaming + mirror
    (3) stella explicit """
    make_comparison_beta0(diff=False)
    make_comparison_beta0(diff=True)

    return

if __name__ == "__main__":
    for_thesis_beta004_benchmark_plot()
    for_thesis_beta0_benchmark_plot()
