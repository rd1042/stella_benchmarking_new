"""From the .out.nc file of a nonlinear simulation, extract phi2(t)
and phi2(kx, ky, tfinal) and save to a .pickle file"""

import sys
sys.path.append("../postprocessing_tools")
from helper_ncdf_new import view_ncdf_variables, extract_data_from_ncdf_with_xarray
import matplotlib.pyplot as plt
import numpy as np
import pickle
import re
from scipy.interpolate import griddata


def postprocess_nonlinear_outnc_sim(outnc_longname):
    """ """

    def make_phi2_kxky_plot(zonal=False):
        """ """
        fig = plt.figure(figsize=[12, 12])

        col1_lhs = 0.1
        row1_bottom = 0.1
        plot_dim_x = 0.75
        plot_dim_y = plot_dim_x
        col1_cb_lhs = 0.9
        cb_width = 0.05

        ax1 = fig.add_axes([col1_lhs, row1_bottom, plot_dim_x, plot_dim_y])
        cbax1 = fig.add_axes([col1_cb_lhs, row1_bottom, cb_width, plot_dim_y])

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

        phi2_contours = ax1.contourf(new_kx, new_ky, logphi2_grid, 20, cmap="viridis")
        fig.colorbar(phi2_contours, cax=cbax1)

        ax1.set_xlabel(r"$k_x$")
        ax1.set_ylabel(r"$k_y$")
        fig.suptitle(sim_longname + "\n" + r"$\|\phi|^2$")

        if zonal:
            save_name = sim_longname + "_phi2_kxky.png"
        else:
            save_name = sim_longname + "_phi2_kxky_nonzonal.png"

        plt.savefig(save_name)
        plt.close()

        return

    [t, kx, ky, phi2, phi2_tkxky] = extract_data_from_ncdf_with_xarray(outnc_longname,
                                    "t", "kx", "ky", "phi2", "phi2_vs_kxky")
    t = np.array(t)
    kx = np.array(kx)
    ky = np.array(ky)
    phi2 = np.array(phi2)
    phi2_tkxky = np.array(phi2_tkxky)

    sim_longname = re.split(".out.nc", outnc_longname)[0]
    pickle_longname = sim_longname + ".summary_pickle"

    myfile = open(pickle_longname, "wb")
    pickle.dump([t, kx, ky, phi2, phi2_tkxky[-1,:,:]], myfile)
    myfile.close()

    ### Make some plots

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(t, phi2_t, c="black", lw=2)
    ax1.set_yscale("log")
    ax1.set_xlabel("t")
    ax1.set_ylabel("phi2")
    ax1.grid(True)
    plt.savefig(sim_longname + "_phi2_t.eps")
    plt.close()

    phi2_kxky = phi2_tkxky[-1,:,:]
    make_phi2_kxky_plot(zonal=True)
    make_phi2_kxky_plot(zonal=False)

    return

def postprocess_folder(folder_longname):
    """ """
    outnc_longname = folder_longname + "/input.out.nc"
    postprocess_nonlinear_outnc_sim(outnc_longname)

    return

if __name__ == "__main__":
    postprocess_nonlinear_outnc_sim("sims/test_em_nonlinear/input_master_like_nonlinear_dt0002_fg_full_krange.out.nc")
    postprocess_nonlinear_outnc_sim("sims/test_em_nonlinear/input_master_like_nonlinear_dt001_fg_full_krange.out.nc")
