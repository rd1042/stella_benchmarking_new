""" """

import sys
sys.path.append("../postprocessing_tools")
from plotting_helper import make_beta_scan_plots, make_comparison_plots, plot_gmvus, plot_gzvs
from helper_ncdf import view_ncdf_variables, extract_data_from_ncdf
from extract_sim_data import find_bpar_phi_ratio, find_apar_phi_ratio
import matplotlib.pyplot as plt
import numpy as np
import glob
import re
import pickle
import xarray as xr
from scipy.interpolate import griddata
# import ast

IMAGE_DIR = "./images/"

leapfrog_nonlinear_outnc = "stella_nonlinear_adiabatic_leapfrog/input.out.nc"
master_nonlinear_outnc = "stella_nonlinear_adiabatic_master/input.out.nc"


def plot_phi2_spectrum(outnc_file):
    """ """

    # [kx, ky, phi2_tkxky] = extract_from_outnc_file(outnc_file, ["kx", "ky", "phi2_vs_kxky"])
    # # Find the latest time at which phi2 has no NaN values
    # is_all_finite = False
    # time_idx = -1
    #
    # while not is_all_finite:
    #     print("time_idx = ", time_idx)
    #     phi2_kxky = phi2_tkxky[time_idx,:,:]
    #     is_all_finite = np.isfinite(phi2_kxky).all()
    #     time_idx -= 1

    data = xr.open_dataset(outnc_file)
    print("data.keys = ", data.keys())
    phi2_tkxky = np.array(data.phi2_vs_kxky)
    phi2_kxky = phi2_tkxky[-1,:,:]
    print("phi2_tkxky.shape = ", phi2_tkxky.shape)
    phi2 = np.array(data.phi2)
    t = np.array(data.t)
    kx = np.array(data.kx)
    ky = np.array(data.ky)
    data.close()
    #sys.exit()
    print("kx, ky = ", kx, ky)
    phi_nonzonal_sum = []
    phi_zonal_sum = []
    phi2_kxky[0,0]  = np.nan
    logphi2_list = []
    kx_list = []
    ky_list = []
    print("kx = ", sorted(set(kx)))
    print("ky = ", sorted(set(ky)))

    for i in range(0, len(kx)):
        for j in range(0, len(ky)):
            kx_list.append(kx[i])
            ky_list.append(ky[j])
            logphi2_list.append(np.log(phi2_kxky[i,j]))

    [kx_min, kx_max] = [min(kx), max(kx)]
    [ky_min, ky_max] = [min(ky), max(ky)]

    new_kx = np.linspace(kx_min, kx_max, 100)
    new_ky = np.linspace(ky_min, ky_max, 100)

    # Turn into a "meshgrid"
    new_kx, new_ky = np.meshgrid(new_kx, new_ky)

    interpolated_data = []  # To store output.

    logphi2_grid = griddata((kx_list, ky_list),
                          np.asarray(logphi2_list), (new_kx, new_ky),
                          method="nearest", fill_value=np.nan)

    fig = plt.figure(figsize=[12, 12])

    col1_lhs = 0.1
    row1_bottom = 0.1
    plot_dim_x = 0.75
    plot_dim_y = plot_dim_x
    col1_cb_lhs = 0.9
    cb_width = 0.05

    ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
    cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])


    phi2_contours = ax1.contourf(new_kx, new_ky, logphi2_grid, 20, cmap="viridis")
    fig.colorbar(phi2_contours, cax=cbax1)

    ax1.set_xlabel(r"$k_x$")
    ax1.set_ylabel(r"$k_y$")
    fig.suptitle(outnc_file + "\n" + r"$\|\phi|^2$")

    plt.show()


    #####################################################
    # Dirty hack to make the same plot as above, but missing the ky=0 row
    # Motivation: It *might* be physically reasonable for most of the energy to
    # be dumped in ky=0 modes (cf. zonal flows)
    #######################################################

    logphi2_list = []
    kx_list = []
    ky_list = []

    for i in range(0, len(kx)):
        for j in range(1, len(ky)): # Adjusted ky limit
            kx_list.append(kx[i])
            ky_list.append(ky[j])
            logphi2_list.append(np.log(phi2_kxky[i,j]))

    [kx_min, kx_max] = [min(kx), max(kx)]
    [ky_min, ky_max] = [min(ky[1:]), max(ky[1:])]
    #print("ky_list = ", ky_list)
    # Define new fprim, tprim coordinates
    new_kx = np.linspace(kx_min, kx_max, 100)
    new_ky = np.linspace(ky_min, ky_max, 100)

    # Turn into a "meshgrid"
    new_kx, new_ky = np.meshgrid(new_kx, new_ky)
    #print("new_ky = ", new_ky)
    interpolated_data = []  # To store output.

    logphi2_grid = griddata((kx_list, ky_list),
                          np.asarray(logphi2_list), (new_kx, new_ky),
                          method="nearest", fill_value=np.nan)


    fig = plt.figure(figsize=[12, 12])

    col1_lhs = 0.1
    row1_bottom = 0.1
    plot_dim_x = 0.75
    plot_dim_y = plot_dim_x
    col1_cb_lhs = 0.9
    cb_width = 0.05

    ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
    cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])


    phi2_contours = ax1.contourf(new_kx, new_ky, logphi2_grid, 20, cmap="viridis")
    fig.colorbar(phi2_contours, cax=cbax1)

    ax1.set_xlabel(r"$k_x$")
    ax1.set_ylabel(r"$k_y$")
    fig.suptitle(outnc_file + "\n" + r"$\log(|\phi|^2)$")

    plt.show()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(t, phi2)
    ax1.grid(True)
    plt.show()


    return


if __name__ == "__main__":
    print("Hello world")
    #view_ncdf_variables(leapfrog_nonlinear_outnc)

    #plot_phi2_spectrum(leapfrog_nonlinear_outnc)
    plot_phi2_spectrum(master_nonlinear_outnc)

    print("Done")
