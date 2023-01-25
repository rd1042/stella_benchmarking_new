""" """

import sys
sys.path.append("../postprocessing_tools")
import matplotlib.scale as scale
import matplotlib.pyplot as plt
import numpy as np
import glob
import re
import pickle
from scipy.interpolate import griddata

def get_data_from_pickle(pickle_longname, fluxes=False):
    """ """

    myfile = open(pickle_longname, "rb")
    if fluxes:
        [t, kx, ky, phi2, phi2_tkxky, parsed_pflx, parsed_vflx, parsed_qflx] = pickle.load(myfile)
    else:
        [t, kx, ky, phi2, phi2_tkxky] = pickle.load(myfile)

    myfile.close()

    sort_idxs = np.argsort(kx)
    kx = kx[sort_idxs]
    print("phi2_tkxky.shape = ", phi2_tkxky.shape)
    phi2_tkxky = phi2_tkxky[:,sort_idxs,:]

    if fluxes:
        parsed_pflx = parsed_pflx[:,:,sort_idxs,:]
        parsed_vflx = parsed_vflx[:,:,sort_idxs,:]
        parsed_qflx = parsed_qflx[:,:,sort_idxs,:]
        return [t, kx, ky, phi2, phi2_tkxky, parsed_pflx, parsed_vflx, parsed_qflx]
    else:
        return [t, kx, ky, phi2, phi2_tkxky]

def make_nl_sim_comparison():
    """ """

    master_cfl025_pickle = "../test_cbc_nonlinear_beta0/sims/stella_nonlinear_2species_master_cflcushion025/input.summary_pickle"
    em_dt001_pickle = "sims/test_em_nonlinear/input_master_like_nonlinear_dt001_fg_full_krange.summary_pickle"
    em_dt0002_pickle = "sims/test_em_nonlinear/input_master_like_nonlinear_dt0002_fg_full_krange.summary_pickle"


    pickle_names = [master_cfl025_pickle, em_dt001_pickle, em_dt0002_pickle]
    labels_list=["master", "em, dt=0.01", "em, dt=0.002"]

    marker_size = 5
    linewidth=2
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    
    flux_t = [True, False, False]

    for idx, pickle in enumerate(pickle_names):
        if flux_t[idx] == True:
            [t, kx, ky, phi2, phi2_tkxky, parsed_pflx, parsed_vflx, parsed_qflx] = get_data_from_pickle(pickle, fluxes=True)
        else:
            [t, kx, ky, phi2, phi2_tkxky] = get_data_from_pickle(pickle, fluxes=False)

        ax1.plot(t, phi2, lw=linewidth, zorder=0, label=labels_list[idx])

    ## Dimensions of phi2_tkxky are [t, kx, ky]
    ## Dimensions of parsed_pflx are [t_idxs, spec, kx, ky]

    # ## imshow makes the plot, but in index-space rather than (kx, ky)
    # # To have meaningful labels, we want to find indexes corresponding to the
    # # desired values of kx, ky
    # kx_label_vals = np.array([-2, -1, 0, 1, 2])
    # ky_label_vals = np.array([0, 1, 2])
    # # idx space goes [0, 1, 2, . . . , nkx(y)-1]
    # # kx(y) space goes [kx(y)_min,  . . . , kx(y)_max]
    # # so kx = ((kx_max - kx_min)/(nkx-1))*kxidx + kx_min
    # # so kxidx = (kx - kx_min) * ((nkx-1)/(kx_max - kx_min))
    # kx_min = np.min(kx) ; kx_max = np.max(kx)
    # ky_min = np.min(ky) ; ky_max = np.max(ky)

    # kx_extent = kx_max - kx_min
    # ky_extent = ky_max - ky_min
    # nkx = len(kx) ; nky = len(ky)
    # xlabel_idx_vals = (kx_label_vals - kx_min) * ((nkx-1)/kx_extent)
    # ylabel_idx_vals = (ky_label_vals - ky_min) * ((nky-1)/ky_extent)
    # # print("xlabel_idx_vals = ", xlabel_idx_vals)
    # # print("ylabel_idx_vals = ", ylabel_idx_vals)
    # # sys.exit()

    ax1.legend(loc="best")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$t$")
    ax1.set_ylabel(r"$\phi$^2")
    plt.tight_layout()
    plt.show()

    return






if __name__ == "__main__":
    make_nl_sim_comparison()
