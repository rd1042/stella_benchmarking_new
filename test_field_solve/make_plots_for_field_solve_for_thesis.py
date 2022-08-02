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


def compare_field_for_thesis(field_name, outnc_longnames, sim_types, sim_labels, normalise=False,
                  title=None):
    """ """
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

    linestyles=((0,(1,0)), (0, (2,2)), (0,(5,1,1,1)), (0,(1,1)),
                (0, (3,2,2,3)), (0, (1,0)) )
    linewidths = [6, 4, 3, 3, 3, 2]

    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    z_orig = np.array(z_list[0])
    field_real = np.array(field_real_list[0])
    field_imag = np.array(field_imag_list[0])
    for sim_idx, sim_label in enumerate(sim_labels_list):
        z = np.array(z_list[sim_idx])
        field_real = np.array(field_real_list[sim_idx])
        field_imag = np.array(field_imag_list[sim_idx])
        if normalise:
            norm_fac_real = max(np.max(abs(field_real)), 1E-12)
            field_real = field_real/norm_fac_real
            norm_fac_imag = max(np.max(abs(field_imag)), 1E-12)
            field_imag = field_imag/norm_fac_imag
        ax1.plot(z/np.pi, field_real, ls=linestyles[sim_idx], lw=linewidths[sim_idx], label=sim_label)
        ax3.plot(z/np.pi, field_imag, ls=linestyles[sim_idx], lw=linewidths[sim_idx])

        if len(z_orig) == len(z):
            if np.max(abs(z_orig - z)) < 1E-4:
                ax2.plot(z/np.pi, ((field_real - field_real_list[0])/field_real_list[0])*100, label="% diff", ls=linestyles[sim_idx], lw=linewidths[sim_idx])
                ax4.plot(z/np.pi, ((field_imag - field_imag_list[0])/field_imag_list[0])*100, ls=linestyles[sim_idx], lw=linewidths[sim_idx])
            else:
                print("z vals don't match")
                print("z_orig - z = ", z_orig - z)
                print("np.max(abs(z_orig - z)) = ", np.max(abs(z_orig - z)))
                print("z_orig = ", np.array(z_orig))
                print("z = ", np.array(z))
        else:
            print("z lengths not equal")

    ax1.legend(loc="best")
    ax3.legend(loc="best")
    ax1.set_ylabel("real(" + field_name + ")")
    ax3.set_ylabel("imag(" + field_name + ")")
    ax2.set_ylabel("real(" + field_name + ") % error")
    ax4.set_ylabel("imag(" + field_name + ") % error")
    ax3.set_xlabel(r"$z/\pi$")
    ax4.set_xlabel(r"$z/\pi$")
    if title is not None:
        fig.suptitle(title)
    plt.tight_layout()
    plt.show()
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
    stella_sim_fapar1_fbpar1= "sims/stella_slab_fapar1_fbpar1/input.out.nc"
    stella_sim_fapar1_fbpar1_lr= "sims/stella_slab_fapar1_fbpar1/input_nvgrid12_nmu6.out.nc"
    stella_sim_fapar1_fbpar1_vlr= "sims/stella_slab_fapar1_fbpar1/input_nvgrid6_nmu3.out.nc"
    stella_sim_fapar0_fbpar1= "sims/stella_slab_fapar0_fbpar1/input.out.nc"
    stella_sim_fapar1_fbpar0= "sims/stella_slab_fapar1_fbpar0/input.out.nc"

    compare_phi_for_thesis([gs2_sim_fapar0_fbpar0, gs2_sim_fapar0_fbpar0_lr, gs2_sim_fapar0_fbpar0_vlr,
                 stella_sim_fapar0_fbpar0, stella_sim_fapar0_fbpar0_lr, stella_sim_fapar0_fbpar0_vlr],
                ["gs2", "gs2", "gs2",
                 "stella", "stella", "stella"],
                ["gs2 (nguass=12, negrid=16)", "gs2 (nguass=6, negrid=8)",
                 "gs2 (nguass=3, negrid=4)",
                 "stella (nvgrid=24, nmu=12)",
                 "stella (nvgrid=12, nmu=6)",
                 "stella (nvgrid=6, nmu=3)",
                 ])

    return

if __name__ == "__main__":
    print("Hello world")
    make_field_solve_test_plots_for_thesis()