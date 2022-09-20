""" """

import sys
sys.path.append("../postprocessing_tools")
import matplotlib.pyplot as plt
import numpy as np
import glob
import re
import pickle
from scipy.interpolate import griddata

IMAGE_DIR = "./images/"

leapfrog_pickle = "sims/stella_nonlinear_2species_leapfrog_nonlinear/input.summary_pickle"
master_pickle = "sims/stella_nonlinear_2species_master_archer2/input.summary_pickle"


def plot_phi2_t(pickle_longnames):
    """ """

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    for pickle_longname in pickle_longnames:
        shortname = re.split("/", pickle_longname)[-2]
        myfile = open(pickle_longname, "rb")
        [t, kx, ky, phi2, phi2_kxky] = pickle.load(myfile)
        myfile.close()
        ax1.plot(t, phi2, label=shortname, lw=2)

    ax1.set_yscale("log")
    ax1.set_xlabel(r"$t$")
    ax1.set_ylabel(r"phi2")
    ax1.legend(loc="best")
    plt.show()

    return

def make_colorplot(fig, ax, cbax, kx, ky, phi2_kxky, zonal=True):
    """ """
    kx_list = []
    ky_list = []
    logphi2_list = []
    if zonal:
        for i in range(0, len(kx)):
            for j in range(0, len(ky)): # Adjusted ky limit
                kx_list.append(kx[i])
                ky_list.append(ky[j])
                logphi2_list.append(np.log(phi2_kxky[i,j]))

        [kx_min, kx_max] = [min(kx), max(kx)]
        [ky_min, ky_max] = [min(ky), max(ky)]
    else:
        for i in range(0, len(kx)):
            for j in range(1, len(ky)):  # Adjusted ky limit
                kx_list.append(kx[i])
                ky_list.append(ky[j])
                logphi2_list.append(np.log(phi2_kxky[i,j]))

        [kx_min, kx_max] = [min(kx), max(kx)]
        [ky_min, ky_max] = [min(ky[1:]), max(ky[1:])]

    new_kx = np.linspace(kx_min, kx_max, 100)
    new_ky = np.linspace(ky_min, ky_max, 100)

    # Turn into a "meshgrid"
    new_kx, new_ky = np.meshgrid(new_kx, new_ky)

    interpolated_data = []  # To store output.

    logphi2_grid = griddata((kx_list, ky_list),
                          np.asarray(logphi2_list), (new_kx, new_ky),
                          method="nearest", fill_value=np.nan)

    phi2_contours = ax.contourf(new_kx, new_ky, logphi2_grid, 20, cmap="viridis")
    fig.colorbar(phi2_contours, cax=cbax)

    return


def plot_properties_master_sim():
    """ """
    pickle_longname = master_pickle
    myfile = open(pickle_longname, "rb")
    [t, kx, ky, phi2, phi2_kxky] = pickle.load(myfile)
    myfile.close()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(t, phi2)
    ax1.scatter(t, phi2, c="red")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$t$")
    ax1.set_ylabel(r"phi2")
    plt.show()


    fig = plt.figure(figsize=[12, 12])

    col1_lhs = 0.1
    row1_bottom = 0.1
    plot_dim_x = 0.75
    plot_dim_y = plot_dim_x
    col1_cb_lhs = 0.9
    cb_width = 0.05

    ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
    cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])

    make_colorplot(fig, ax1, cbax1, kx, ky, phi2_kxky, zonal=True)

    ax1.set_xlabel(r"$k_x$")
    ax1.set_ylabel(r"$k_y$")
    plt.show()


    return

if __name__ == "__main__":
    print("Hello world")
    plot_properties_master_sim()
