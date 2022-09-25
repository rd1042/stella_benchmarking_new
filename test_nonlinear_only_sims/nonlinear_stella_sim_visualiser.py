""" """

import sys
sys.path.append("../postprocessing_tools")
from plotting_helper import make_comparison_plots, plot_gmvus, plot_gzvs
from plotting_helper import make_comparison_plots_leapfrog_poster
from helper_ncdf import view_ncdf_variables, extract_data_from_ncdf
from extract_sim_data import get_omega_data
import matplotlib.pyplot as plt
import numpy as np
from numpy import fft
import glob

def make_phi_kxky_modes_pics(outnc_longname, log=False):
    """ """
    # Get phi(kx, ky, z, t)
    [t, kx, ky, z, phi_vs_t] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi_vs_t")
    # print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
    # print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))
    # print("t = ", t)
    # print("len(t) = ", len(t))
    # sys.exit()
    nz_per_mode = len(z)
    z_idx = int(nz_per_mode/2)  # The z idx of z=0
    if z[z_idx] != 0:
        print("Error! z[z_idx] = ", z[z_idx])
        sys.exit()

    # z_idx = 1 # Old
    counter = 0
    phi_t_ky_kx = phi_vs_t[:, 0, z_idx, :, :]
    for ky_idx in range(0, len(ky)):
        for kx_idx in range(0, len(kx)):
            phi_t = phi_t_ky_kx[:, kx_idx, ky_idx]

            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            # ax2 = fig.add_subplot(312, sharex=ax1)
            # ax3 = fig.add_subplot(313, sharex=ax1)
            ax1.plot(t, phi_t.real, label="real(phi)")
            ax1.plot(t, phi_t.imag, label="im(phi)")
            ax1.plot(t, abs(phi_t), label="abs(phi)")
            if log:
                ax1.set_yscale("log")
            ax1.set_xlabel(r"$t$")
            ax1.set_ylabel(r"$\phi$")
            ax1.legend(loc="best")
            fig.suptitle("kx={:.3f}, ky={:.3f}".format(kx[kx_idx], ky[ky_idx]))
            save_name="images/phi_t_{:02d}.png".format(counter)
            counter+=1
            plt.savefig(save_name)
            plt.close()

    return

def make_phi2_kxky_modes_pics_single_mode(outnc_longname):
    """ """
    # Get phi(kx, ky, z, t)
    [t, kx, ky, z, phi_vs_t] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi_vs_t")
    # print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
    # print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))
    # print("kx = ", kx)
    # print("ky = ", ky)
    # print("t = ", t)
    # print("len(t) = ", len(t))
    # sys.exit()
    nz_per_mode = len(z)
    z_idx = int(nz_per_mode/2)  # The z idx of z=0
    if z[z_idx] != 0:
        print("Error! z[z_idx] = ", z[z_idx])
        sys.exit()

    # z_idx = 1 # Old
    counter = 0
    phi_t_ky_kx = phi_vs_t[:, 0, z_idx, :, :]
    for ky_idx in [2]: #range(0, len(ky)):
        for kx_idx in [2]: #range(0, len(kx)):
            phi_t = phi_t_ky_kx[:, kx_idx, ky_idx]

            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            # ax2 = fig.add_subplot(312, sharex=ax1)
            # ax3 = fig.add_subplot(313, sharex=ax1)
            ax1.plot(t, phi_t.real, label="real(phi)")
            ax1.plot(t, phi_t.imag, label="im(phi)")
            ax1.plot(t, abs(phi_t), label="abs(phi)")
            ax1.scatter(t, phi_t.real, marker="x", c="black")
            ax1.scatter(t, phi_t.imag, marker="x", c="black")
            ax1.scatter(t, abs(phi_t), marker="x", c="black")

            ax1.set_xlabel(r"$t$")
            ax1.set_ylabel(r"$\phi$")
            ax1.legend(loc="best")
            fig.suptitle("kx={:.3f}, ky={:.3f}".format(kx[kx_idx], ky[ky_idx]))
            save_name="images/phi_t_{:02d}.png".format(counter)
            counter+=1
            plt.show()

    return

def make_phi2_kxky_modes_pics_single_mode_multiple_sims(outnc_longnames, sim_labels, linestyles):
    """ """
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    for sim_idx, outnc_longname in enumerate(outnc_longnames):
        print("outnc_longname = ", outnc_longname)
        # Get phi(kx, ky, z, t)
        [t, kx, ky, z, phi_vs_t] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi_vs_t")
        # print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
        # print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))
        # print("kx = ", kx)
        # print("ky = ", ky)
        # print("t = ", t)
        # print("len(t) = ", len(t))
        # sys.exit()
        nz_per_mode = len(z)
        z_idx = int(nz_per_mode/2)  # The z idx of z=0
        if z[z_idx] != 0:
            print("Error! z[z_idx] = ", z[z_idx])
            sys.exit()

        # z_idx = 1 # Old
        counter = 0
        phi_t_ky_kx = phi_vs_t[:, 0, z_idx, :, :]
        for ky_idx in [2]: #range(0, len(ky)):
            for kx_idx in [2]: #range(0, len(kx)):
                phi_t = phi_t_ky_kx[:, kx_idx, ky_idx]


                # ax2 = fig.add_subplot(312, sharex=ax1)
                # ax3 = fig.add_subplot(313, sharex=ax1)
                ax1.plot(t, abs(phi_t), label=sim_labels[sim_idx], ls=linestyles[sim_idx])
                # ax1.scatter(t, phi_t.real, marker="x", c="black")
                # ax1.scatter(t, phi_t.imag, marker="x", c="black")
                #ax1.scatter(t, abs(phi_t), marker="x", c="black")

    ax1.set_xlabel(r"$t$")
    ax1.set_ylabel(r"$\phi$")
    ax1.legend(loc="best")
    fig.suptitle("kx={:.3f}, ky={:.3f}".format(kx[kx_idx], ky[ky_idx]))
    save_name="images/phi_t_{:02d}.png".format(counter)
    counter+=1
    plt.show()

    return

def make_phi_kxky_modes_pics_multiple_sims(outnc_longnames, sim_labels, linestyles):
    """ """
    t_list = []
    kx_list = []
    ky_list = []
    phi_t_ky_kx_list = []

    for sim_idx, outnc_longname in enumerate(outnc_longnames):
        print("outnc_longname = ", outnc_longname)
        # Get phi(kx, ky, z, t)
        [t, kx, ky, z, phi_vs_t] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi_vs_t")
        t_list.append(t)
        kx_list.append(kx)
        ky_list.append(ky)
        nz_per_mode = len(z)
        z_idx = int(nz_per_mode/2)  # The z idx of z=0
        if z[z_idx] != 0:
            print("Error! z[z_idx] = ", z[z_idx])
            sys.exit()
        phi_t_ky_kx_list.append(phi_vs_t[:, 0, z_idx, :, :])
        #print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
        # print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))
        print("kx = ", kx)
        print("ky = ", ky)
        # print("t = ", t)
        # print("len(t) = ", len(t))
        # sys.exit()
    ky_vals = sorted(set(np.concatenate(ky_list)))
    kx_vals = sorted(set(np.concatenate(kx_list)))
    print("ky_vals = ", ky_vals)
    print("kx_vals = ", kx_vals)
    # sys.exit()
    counter = 0
    for ky_val in ky_vals :
        for kx_val in kx_vals:
            fig = plt.figure(figsize=[14,10])
            ax1 = fig.add_subplot(111)
            for sim_idx in range(0, len(outnc_longnames)):
                ky = ky_list[sim_idx]; kx = kx_list[sim_idx]
                ky_idx = np.argmin(abs(ky-ky_val))
                kx_idx = np.argmin(abs(kx-kx_val))
                if (abs(ky[ky_idx] - ky_val) < 0.001) and (abs(kx[kx_idx] - kx_val) < 0.001):
                    t = t_list[sim_idx]
                    phi_t_ky_kx = phi_t_ky_kx_list[sim_idx]
                    phi_t = phi_t_ky_kx[:, kx_idx, ky_idx]


                    # ax2 = fig.add_subplot(312, sharex=ax1)
                    # ax3 = fig.add_subplot(313, sharex=ax1)
                    ax1.plot(t, abs(phi_t), label=sim_labels[sim_idx], ls=linestyles[sim_idx], lw=3)
                    # ax1.scatter(t, phi_t.real, marker="x", c="black")
                    # ax1.scatter(t, phi_t.imag, marker="x", c="black")
                    #ax1.scatter(t, abs(phi_t), marker="x", c="black")

            ax1.set_yscale("log")
            ax1.set_xlabel(r"$t$")
            ax1.set_ylabel(r"$\phi$")
            ax1.legend(loc="best")
            fig.suptitle("kx={:.3f}, ky={:.3f}".format(kx_val, ky_val))
            save_name="images/phi_t_{:02d}.png".format(counter)
            plt.savefig(save_name)
            plt.close()
            counter+=1

    return

def make_phi_kxky_modes_pics_multiple_sims_all_z(outnc_longnames, sim_labels, cols, transparency_val=0.5):
    """As above, but plotting phi for all z values """
    t_list = []
    kx_list = []
    ky_list = []
    phi_t_z_ky_kx_list = []

    for sim_idx, outnc_longname in enumerate(outnc_longnames):
        print("outnc_longname = ", outnc_longname)
        # Get phi(kx, ky, z, t)
        [t, kx, ky, z, phi_vs_t] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi_vs_t")
        t_list.append(t)
        kx_list.append(kx)
        ky_list.append(ky)
        nz_per_mode = len(z)
        phi_t_z_ky_kx_list.append(phi_vs_t[:, 0, :, :, :])
        #print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
        # print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))
        print("kx = ", kx)
        print("ky = ", ky)
        # print("t = ", t)
        # print("len(t) = ", len(t))
        # sys.exit()
    ky_vals = sorted(set(np.concatenate(ky_list)))
    kx_vals = sorted(set(np.concatenate(kx_list)))
    print("ky_vals = ", ky_vals)
    print("kx_vals = ", kx_vals)
    sys.exit()
    counter = 0
    for ky_val in ky_vals :
        for kx_val in kx_vals:
            fig = plt.figure(figsize=[14,10])
            ax1 = fig.add_subplot(111)
            for sim_idx in range(0, len(outnc_longnames)):
                ky = ky_list[sim_idx]; kx = kx_list[sim_idx]
                ky_idx = np.argmin(abs(ky-ky_val))
                kx_idx = np.argmin(abs(kx-kx_val))
                if (abs(ky[ky_idx] - ky_val) < 0.001) and (abs(kx[kx_idx] - kx_val) < 0.001):
                    t = t_list[sim_idx]
                    phi_t_z_ky_kx = phi_t_z_ky_kx_list[sim_idx]
                    phi_t_z = phi_t_z_ky_kx[:, :, kx_idx, ky_idx]


                    # ax2 = fig.add_subplot(312, sharex=ax1)
                    # ax3 = fig.add_subplot(313, sharex=ax1)
                    ax1.plot(t, abs(phi_t_z[:,0]), label=sim_labels[sim_idx], c=cols[sim_idx], lw=3, alpha=transparency_val)
                    for zidx in range(len(phi_t_z[0,:])):
                        ax1.plot(t, abs(phi_t_z[:,zidx]), c=cols[sim_idx], lw=3, alpha=transparency_val)
                    # ax1.scatter(t, phi_t.real, marker="x", c="black")
                    # ax1.scatter(t, phi_t.imag, marker="x", c="black")
                    #ax1.scatter(t, abs(phi_t), marker="x", c="black")

            ax1.set_yscale("log")
            ax1.set_xlabel(r"$t$")
            ax1.set_ylabel(r"$\phi$")
            ax1.legend(loc="best")
            fig.suptitle("kx={:.3f}, ky={:.3f}".format(kx_val, ky_val))
            save_name="images/phi_t_{:02d}.png".format(counter)
            plt.savefig(save_name)
            plt.close()
            counter+=1

    return

def make_phi2_kxky_modes_pics_multiple_sims(outnc_longnames, sim_labels, linestyles):
    """ """
    t_list = []
    kx_list = []
    ky_list = []
    phi2_kxky_list = []

    for sim_idx, outnc_longname in enumerate(outnc_longnames):
        print("outnc_longname = ", outnc_longname)
        # Get phi(kx, ky, z, t)
        [t, kx, ky, z, phi2_kxky] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi2_vs_kxky")
        t_list.append(t)
        kx_list.append(kx)
        ky_list.append(ky)
        nz_per_mode = len(z)

        phi2_kxky_list.append(phi2_kxky)
        #print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
        # print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))
        # print("kx = ", kx)
        # print("ky = ", ky)
        # print("t = ", t)
        # print("len(t) = ", len(t))
        # sys.exit()
    ky_vals = sorted(set(np.concatenate(ky_list)))
    kx_vals = sorted(set(np.concatenate(kx_list)))
    print("ky_vals = ", ky_vals)
    print("kx_vals = ", kx_vals)
    # sys.exit()
    counter = 0
    for ky_val in ky_vals :
        for kx_val in kx_vals:
            fig = plt.figure(figsize=[14,10])
            ax1 = fig.add_subplot(111)
            for sim_idx in range(0, len(outnc_longnames)):
                ky = ky_list[sim_idx]; kx = kx_list[sim_idx]
                ky_idx = np.argmin(abs(ky-ky_val))
                kx_idx = np.argmin(abs(kx-kx_val))
                if (abs(ky[ky_idx] - ky_val) < 0.001) and (abs(kx[kx_idx] - kx_val) < 0.001):
                    t = t_list[sim_idx]
                    phi2_kxky = phi2_kxky_list[sim_idx]
                    phi2_t = phi2_kxky[:, kx_idx, ky_idx]


                    # ax2 = fig.add_subplot(312, sharex=ax1)
                    # ax3 = fig.add_subplot(313, sharex=ax1)
                    ax1.plot(t, phi2_t, label=sim_labels[sim_idx], ls=linestyles[sim_idx], lw=3)
                    # ax1.scatter(t, phi_t.real, marker="x", c="black")
                    # ax1.scatter(t, phi_t.imag, marker="x", c="black")
                    #ax1.scatter(t, abs(phi_t), marker="x", c="black")

            ax1.set_yscale("log")
            ax1.set_xlabel(r"$t$")
            ax1.set_ylabel(r"phi2")
            ax1.legend(loc="best")
            fig.suptitle("kx={:.3f}, ky={:.3f}".format(kx_val, ky_val))
            save_name="images/phi_t_{:02d}.png".format(counter)
            plt.savefig(save_name)
            plt.close()
            counter+=1

    return

def plot_phiz_for_zonal_mode(outnc_longname):
    """For zonal modes, plot phi(z); is it periodic? """

    # view_ncdf_variables(outnc_longname)
    # sys.exit()
    [t, kx, ky, z, phi_vs_t] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi_vs_t")
    #print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky

    print("max(abs(phi_vs_t[-1,0,-1,:,0] - phi_vs_t[-1,0,0,:,0])) = ", np.max(abs(phi_vs_t[-1,0,-1,:,0] - phi_vs_t[-1,0,0,:,0])))
    print("max(abs(phi_vs_t[-1,0,-1,:,:] - phi_vs_t[-1,0,0,:,:])) = ", np.max(abs(phi_vs_t[-1,0,-1,:,:] - phi_vs_t[-1,0,0,:,:])))
    #sys.exit()
    z_plus_max_z = z + np.max(z) - np.min(z)
    duplicated_z = np.concatenate((z, z_plus_max_z))
    for ky_idx in [0,1, 10, 20, 30]:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        for kx_idx in range(0, len(kx)):
            phi_z_norm = phi_vs_t[-1, 0, :, kx_idx, ky_idx]/np.max(abs(phi_vs_t[-1, 0, :, kx_idx, ky_idx]))
            duplicated_phiz = np.concatenate((phi_z_norm, phi_z_norm))
            ax1.plot(duplicated_z/np.pi, abs(duplicated_phiz))
            ax1.scatter(duplicated_z/np.pi, abs(duplicated_phiz), marker="x")
        ax1.plot([np.max(z)/np.pi, np.max(z)/np.pi], [0, 1.1], c="black", ls="--", lw=2, alpha=0.5)
        #ax1.set_yscale("log")
        ax1.set_xlabel(r"$z/\pi$")
        ax1.set_ylabel(r"$\phi$")
        plt.show()
    return

def make_phi_ky_modes_pics_multiple_sims(outnc_longnames, sim_labels, sim_cols,
                    min_phi=1E-5, transparency_val=0.1):
    """ """
    t_list = []
    kx_list = []
    ky_list = []
    phi_t_ky_kx_list = []

    for sim_idx, outnc_longname in enumerate(outnc_longnames):
        print("outnc_longname = ", outnc_longname)
        # Get phi(kx, ky, z, t)
        [t, kx, ky, z, phi_vs_t] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi_vs_t")
        nz_per_mode = len(z)
        z_idx = int(nz_per_mode/2)  # The z idx of z=0
        t_list.append(t)
        kx_list.append(kx)
        ky_list.append(ky)
        nz_per_mode = len(z)

        phi_t_ky_kx_list.append(phi_vs_t[:, 0, z_idx, :, :])
        #print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
        # print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))
        # print("kx = ", kx)
        # print("ky = ", ky)
        # print("t = ", t)
        # print("len(t) = ", len(t))
        # sys.exit()
    ky_vals = sorted(set(np.concatenate(ky_list)))
    kx_vals = sorted(set(np.concatenate(kx_list)))
    print("ky_vals = ", ky_vals)
    print("kx_vals = ", kx_vals)
    # sys.exit()
    counter = 0
    for ky_val in ky_vals :
        fig = plt.figure(figsize=[16,14])
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        for sim_idx in range(0, len(outnc_longnames)):
            ky = ky_list[sim_idx]; kx = kx_list[sim_idx]
            ky_idx = np.argmin(abs(ky-ky_val))
            max_phi = 0
            if (abs(ky[ky_idx] - ky_val) < 0.001):
                t = t_list[sim_idx]
                phi_kxky = phi_t_ky_kx_list[sim_idx]
                phi_kx0_t = phi_kxky[:, 0, ky_idx]
                ax1.plot(t, abs(phi_kx0_t), c=sim_cols[sim_idx], lw=3, alpha=transparency_val, label=sim_labels[sim_idx])
                kx = kx_list[sim_idx]
                sorted_kx_idxs = np.argsort(kx)
                ax2.plot(kx[sorted_kx_idxs], abs(phi_kxky[-1,sorted_kx_idxs,ky_idx]), c=sim_cols[sim_idx], lw=3, label=sim_labels[sim_idx])
                max_phi = max(max_phi, np.max(abs(phi_kx0_t)))
                for kx_idx in range(1, len(kx_list[sim_idx])):
                    phi_t = phi_kxky[:, kx_idx, ky_idx]
                    ax1.plot(t, abs(phi_t), c=sim_cols[sim_idx], lw=3, alpha=transparency_val)
                    max_phi = max(max_phi, np.max(abs(phi_t)))

                    # ax2 = fig.add_subplot(312, sharex=ax1)
                    # ax3 = fig.add_subplot(313, sharex=ax1)
                    # ax1.scatter(t, phi_t.real, marker="x", c="black")
                    # ax1.scatter(t, phi_t.imag, marker="x", c="black")
                    #ax1.scatter(t, abs(phi_t), marker="x", c="black")

        ax1.set_yscale("log")
        ax2.set_yscale("log")
        ax1.set_xlabel(r"$t$")
        ax2.set_xlabel(r"$k_x$")
        ax1.set_ylabel(r"abs(phi)")
        ax2.set_ylabel(r"abs(phi(t=tfinal))")
        ax1.legend(loc="best")
        ax2.legend(loc="best")
        ax1.set_ylim(min_phi, max_phi)
        fig.suptitle("ky={:.3f}".format(ky_val))
        save_name="images/phi_t_{:02d}.png".format(counter)
        plt.savefig(save_name)
        plt.close()
        counter+=1

    return

def make_phi2_ky_modes_pics_multiple_sims(outnc_longnames, sim_labels, sim_cols):
    """ """
    t_list = []
    kx_list = []
    ky_list = []
    phi2_kxky_list = []

    for sim_idx, outnc_longname in enumerate(outnc_longnames):
        print("outnc_longname = ", outnc_longname)
        # Get phi(kx, ky, z, t)
        [t, kx, ky, z, phi2_kxky] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi2_vs_kxky")
        t_list.append(t)
        kx_list.append(kx)
        ky_list.append(ky)
        nz_per_mode = len(z)

        phi2_kxky_list.append(phi2_kxky)
        #print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
        # print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))
        # print("kx = ", kx)
        # print("ky = ", ky)
        # print("t = ", t)
        # print("len(t) = ", len(t))
        # sys.exit()
    ky_vals = sorted(set(np.concatenate(ky_list)))
    kx_vals = sorted(set(np.concatenate(kx_list)))
    print("ky_vals = ", ky_vals)
    print("kx_vals = ", kx_vals)
    # sys.exit()
    counter = 0
    transparency_val = 0.1
    for ky_val in ky_vals :
        fig = plt.figure(figsize=[16,14])
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        for sim_idx in range(0, len(outnc_longnames)):
            ky = ky_list[sim_idx]; kx = kx_list[sim_idx]
            ky_idx = np.argmin(abs(ky-ky_val))
            max_phi2 = 0
            if (abs(ky[ky_idx] - ky_val) < 0.001):
                t = t_list[sim_idx]
                phi2_kxky = phi2_kxky_list[sim_idx]
                phi2_t = phi2_kxky[:, 0, ky_idx]
                ax1.plot(t, phi2_t, c=sim_cols[sim_idx], lw=3, alpha=transparency_val, label=sim_labels[sim_idx])
                kx = kx_list[sim_idx]
                sorted_kx_idxs = np.argsort(kx)
                ax2.plot(kx[sorted_kx_idxs], phi2_kxky[-1,sorted_kx_idxs,ky_idx], c=sim_cols[sim_idx], lw=3, label=sim_labels[sim_idx])
                max_phi2 = max(max_phi2, np.max(phi2_t))
                for kx_idx in range(1, len(kx_list[sim_idx])):
                    phi2_t = phi2_kxky[:, kx_idx, ky_idx]
                    ax1.plot(t, phi2_t, c=sim_cols[sim_idx], lw=3, alpha=transparency_val)
                    max_phi2 = max(max_phi2, np.max(phi2_t))

                    # ax2 = fig.add_subplot(312, sharex=ax1)
                    # ax3 = fig.add_subplot(313, sharex=ax1)
                    # ax1.scatter(t, phi_t.real, marker="x", c="black")
                    # ax1.scatter(t, phi_t.imag, marker="x", c="black")
                    #ax1.scatter(t, abs(phi_t), marker="x", c="black")

        ax1.set_yscale("log")
        ax2.set_yscale("log")
        ax1.set_xlabel(r"$t$")
        ax2.set_xlabel(r"$k_x$")
        ax1.set_ylabel(r"phi2")
        ax2.set_ylabel(r"phi2(t=tfinal)")
        ax1.legend(loc="best")
        ax2.legend(loc="best")
        ax1.set_ylim(1E-10, max_phi2)
        fig.suptitle("ky={:.3f}".format(ky_val))
        save_name="images/phi_t_{:02d}.png".format(counter)
        plt.savefig(save_name)
        plt.close()
        counter+=1

    return

def make_phi2_ky_modes_pics_multiple_sims_for_talk(outnc_longnames, sim_labels, sim_cols):
    """ """
    t_list = []
    kx_list = []
    ky_list = []
    phi2_kxky_list = []

    for sim_idx, outnc_longname in enumerate(outnc_longnames):
        print("outnc_longname = ", outnc_longname)
        # Get phi(kx, ky, z, t)
        [t, kx, ky, z, phi2_kxky] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi2_vs_kxky")
        t_list.append(t)
        kx_list.append(kx)
        ky_list.append(ky)
        nz_per_mode = len(z)

        phi2_kxky_list.append(phi2_kxky)
        #print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
        # print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))
        # print("kx = ", kx)
        # print("ky = ", ky)
        # print("t = ", t)
        # print("len(t) = ", len(t))
        # sys.exit()
    ky_vals = sorted(set(np.concatenate(ky_list)))
    kx_vals = sorted(set(np.concatenate(kx_list)))
    print("len(kx), len(ky) = ", len(kx), len(ky))
    print("ky_vals = ", ky_vals)
    print("kx_vals = ", kx_vals)

    counter = 0
    transparency_val = 0.1
    tfinal =   400
    for ky_val in ky_vals :
        fig = plt.figure(figsize=[16,14])
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        for sim_idx in range(0, len(outnc_longnames)):
            ky = ky_list[sim_idx]; kx = kx_list[sim_idx]
            ky_idx = np.argmin(abs(ky-ky_val))
            max_phi2 = 0
            if (abs(ky[ky_idx] - ky_val) < 0.001):
                t = t_list[sim_idx]
                phi2_kxky = phi2_kxky_list[sim_idx]
                phi2_t = phi2_kxky[:, 0, ky_idx]
                ax1.plot([450, 460], [1e-2, 1e-2], c=sim_cols[sim_idx], lw=3, label=sim_labels[sim_idx])
                kx = kx_list[sim_idx]
                sorted_kx_idxs = np.argsort(kx)
                ax2.plot(kx[sorted_kx_idxs], phi2_kxky[-1,sorted_kx_idxs,ky_idx], c=sim_cols[sim_idx], lw=3, label=sim_labels[sim_idx])
                max_phi2 = max(max_phi2, np.max(phi2_t))
                for kx_idx in range(0, len(kx_list[sim_idx])):
                    phi2_t = phi2_kxky[:, kx_idx, ky_idx]
                    ax1.plot(t, phi2_t, c=sim_cols[sim_idx], lw=3, alpha=transparency_val)
                    max_phi2 = max(max_phi2, np.max(phi2_t))

                    # ax2 = fig.add_subplot(312, sharex=ax1)
                    # ax3 = fig.add_subplot(313, sharex=ax1)
                    # ax1.scatter(t, phi_t.real, marker="x", c="black")
                    # ax1.scatter(t, phi_t.imag, marker="x", c="black")
                    #ax1.scatter(t, abs(phi_t), marker="x", c="black")

        ax1.set_yscale("log")
        ax2.set_yscale("log")
        ax1.set_xlim(0,tfinal)
        ax1.set_ylim(1e-9, 1)
        ax1.set_xlabel(r"$t (r_{LCFS}/v_{th,i})$", fontsize=35)
        ax2.set_xlabel(r"$k_x \rho_i$", fontsize=35)
        ax1.set_ylabel(r"$\vert \phi_{k_x} \vert^2$", fontsize=35)
        ax2.set_ylabel(r"$\vert \phi_{k_x} (t=t_{final})\vert^2$", fontsize=35)
        ax1.legend(loc="best", fontsize=25)
        #ax2.legend(loc="best", fontsize=25)
        ax1.set_ylim(1E-10, max_phi2)
        for ax in [ax1, ax2]:
            ax.tick_params(labelsize=25)
        #fig.suptitle("ky={:.3f}".format(ky_val))
        plt.tight_layout()
        save_name="images/phi_t_{:02d}.png".format(counter)
        plt.savefig(save_name)
        plt.close()
        counter+=1
        # sys.exit()

    return

def make_phi2_kxky_modes_pics_for_each_z(outnc_longname):
    """ """
    # Get phi(kx, ky, z, t)
    [t, kx, ky, z, phi_vs_t] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi_vs_t")
    # print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
    # print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))
    nz_per_mode = len(z)

    # z_idx = 1 # Old
    counter = 0
    for ky_idx in range(0, len(ky)):
        for kx_idx in range(0, len(kx)):

            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            # ax2 = fig.add_subplot(312, sharex=ax1)
            # ax3 = fig.add_subplot(313, sharex=ax1)
            for z_idx in range(0, nz_per_mode):
                ax1.plot(t, abs(phi_vs_t[:, 0, z_idx, kx_idx,ky_idx]))

            ax1.set_xlabel(r"$t$")
            ax1.set_ylabel(r"$\phi$")
            fig.suptitle("kx={:.3f}, ky={:.3f}".format(kx[kx_idx], ky[ky_idx]))
            save_name="images/phi_t_{:02d}.png".format(counter)
            counter+=1
            plt.savefig(save_name)
            plt.close()
            # plt.show()

    return

def examine_initialisation():
    """ """
    outnc_rk3_nz4 = "looking_at_initial_conditions/rk3_nonlinear_only_nz4.out.nc"
    outnc_nisl_nz4 = "looking_at_initial_conditions/nisl_nonlinear_only_nz4.out.nc"
    outnc_nisl_nz2 = "looking_at_initial_conditions/nisl_nonlinear_only_nz2.out.nc"

    # Get phi(kx, ky, z, t)
    [t_nisl, kx_nisl, ky_nisl, z_nisl, phi_vs_t_nisl] = extract_data_from_ncdf(outnc_nisl_nz4, "t", 'kx', 'ky', "zed", "phi_vs_t")
    [t_rk3, kx_rk3, ky_rk3, z_rk3, phi_vs_t_rk3] = extract_data_from_ncdf(outnc_rk3_nz4, "t", 'kx', 'ky', "zed", "phi_vs_t")
    [t_nz2, kx_nz2, ky_nz2, z_nz2, phi_vs_t_nz2] = extract_data_from_ncdf(outnc_nisl_nz2, "t", 'kx', 'ky', "zed", "phi_vs_t")
    # print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
    phi_nisl = phi_vs_t_nisl[0,0,:,:,:]
    phi_rk3 = phi_vs_t_rk3[0,0,:,:,:]
    # phi_diff = phi_nisl - phi_rk3
    # print("phi_diff = ", phi_diff)
    # print("kx_nisl[1], kx_nisl[-1] = ", kx_nisl[1], kx_nisl[-1])
    print("kx_nisl = ", kx_nisl)
    print("ky_nisl = ", ky_nisl)
    sys.exit()
    print("z, nzed=4 =  ", z_nisl)
    print("z, nzed=2 = ", z_nz2)
    print("phi_vs_t_nz2[0,0,:,1,1] = ", phi_vs_t_nz2[0,0,:,1,1])
    print("phi_vs_t_nz2[0,0,:,-1,1] = ", phi_vs_t_nz2[0,0,:,-1,1])
    # print("phi_nisl[:,1,1] = ", phi_nisl[:,1,1])
    # print("phi_nisl[:,-1,1] = ", phi_nisl[:,-1,1])
    # print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))

def ifft_phi(kx, ky, phi_kxky, extra_upsample_fac=1):
    """ """
    # Before FTing, we need to increase the size of phi_kxky; this is because
    # ifft2 is expecting both kx and ky to include negative values i.e.
    # expects kx = (0, dkx, 2*dkx, . . ., kxmax, -kxmax, -kxmax+dkx, . . . , -dkx)
    # and     ky = (0, dky, 2*dky, . . ., kymax, -kymax, -kymax+dky, . . . , -dky)
    # whereas currently ky has only positive values. This is becasue g(x,y) is
    # purely real, and so the entries for the -ve values of ky are just the
    # complex conjugate of the +ve values, so can be ignored.
    # Question: How to do this correctly to ensure g(x,y) purely real? Have
    # tried padding phi_kxky with the complex conjugate of phi_kxky, but
    # this seems not to work. stella performs the transform in a rather curious
    # way, so could look at stella's routine and adopt it for ourselves?
    # stella seems to:
    # (1) "swap" kx,ky. This consists of changing the shape from (naky, nakx)
    # to (2*naky-1, ikx_max); so going from
    #    kx=(0, . . ., kxmax, -kxmax, . . . , -dkx) , ky=(0, . . ., kymax)
    # to
    #    kx=(0, . . ., kxmax) , ky=(0, . . ., kymax, . . . , -dky)
    # The -ve ky values for a particular kx look to be the complex conjugate
    # of the -ve kx vals for that kx (with kx=0 being treated specially.)
    # (2) Perform a complex-to-complex transformation in y. This is straightforward,
    # except for the padding with zeros to avoid de-aliasing.
    # (3) Perform a complex-to-real transformation in x. It also looks like this
    # is just padding with zeros to avoid de-aliasing; the size of the transformed
    # array is nx/2+1 and the number of non-zero entries is (nakx/2+1). HOWEVER,
    # the output of the tranformation is size nx; so this DOESN'T look like an
    # ordinary FT which just throws away the complex part.

    # First, do the swap. Copied from swap_kxky_complex in ktgrids.f90
    nkx_inc_negative = len(kx) ; nky_no_negative = len(ky)
    nkx_no_negative = int((nkx_inc_negative + 1)/2)

    ### If extra padding is required, add it in here
    if extra_upsample_fac > 1:

        phi_kxky_upsampled = np.zeros((int(2*extra_upsample_fac*(nkx_no_negative-1)+1), int(extra_upsample_fac*(nky_no_negative-1)+1)), dtype="complex")
        # Populate +ve kx, +ky vals
        phi_kxky_upsampled[:nkx_no_negative,:nky_no_negative] = phi_kxky[:nkx_no_negative,:]
        # Populate -ve kx, +ky vals
        phi_kxky_upsampled[-nkx_no_negative+1:,:nky_no_negative] = phi_kxky[-nkx_no_negative+1:,:]
        # print("Done!")
        # print("len(phi_kxky[:nkx_no_negative,0]) = ", len(phi_kxky[:nkx_no_negative,0]))
        # print("len(pphi_kxky[-nkx_no_negative+1:,0]) = ", len(phi_kxky[-nkx_no_negative+1:,0]))
        # print("phi_kxky[2,:] = ", phi_kxky[2,:])
        # print("phi_kxky_upsampled[2,:] = ", phi_kxky_upsampled[2,:])
        # print("phi_kxky[:,2] = ", phi_kxky[:,2])
        # print("phi_kxky_upsampled[:,2] = ", phi_kxky_upsampled[:,2])
        # sys.exit()
        # print("")
        ## Update the dimensions
        nkx_inc_negative = int(2*extra_upsample_fac*(nkx_no_negative-1)+1)
        nky_no_negative = int(extra_upsample_fac*(nky_no_negative-1)+1)
        phi_kxky = phi_kxky_upsampled

    # First, do the swap. Copied from swap_kxky_complex in ktgrids.f90
    nkx_no_negative = int((nkx_inc_negative + 1)/2)
    nky_inc_negative = int(nky_no_negative*2 - 1)
    # print("nkx_no_negative, nkx_inc_negative, nky_no_negative, nky_inc_negative = ",
    #         nkx_no_negative, nkx_inc_negative, nky_no_negative, nky_inc_negative)
    # NB phi_kxky is shape(nkx_inc_negative, nky_no_negative)
    phi_kxky_swap = np.zeros((nkx_no_negative, nky_inc_negative), dtype="complex")
    # +ve kx, ky entries are the same
    phi_kxky_swap[:nkx_no_negative, :nky_no_negative] = phi_kxky[:nkx_no_negative, :nky_no_negative]
    # Treat kx=0 specially
    for iky in range(nky_no_negative,nky_inc_negative):
        # Want the ky idx corresponding to the +ve ky value of ky(iky)
        # e.g. for iky=nky_no_negative, ky=-kymax, so want the ky idx corresponding
        # to +kymax (in this example, ikyneg=nky_no_negative-1)
        ikyneg = nky_inc_negative - iky
        #print("iky, ikyneg = ", iky, ikyneg)
        phi_kxky_swap[0,iky] = np.conj(phi_kxky[0,ikyneg])

    # Now map the (-ve kx, +ve ky) values to (+ve kx, -ve ky) values
    for ikx in range(1, nkx_no_negative):
        # Get the kx idx corresponding to the -ve kx value of kx(ikx)
        # e.g. for ikx=nkx_no_negative-1 , kx=kxmax so want the kx idx
        # corresponding to -kxmax (in this case, nkx_no_negative)
        ikxneg = nkx_inc_negative - ikx
        for iky in range(nky_no_negative,nky_inc_negative):
            ikyneg = nky_inc_negative - iky
            phi_kxky_swap[ikx,iky] = np.conj(phi_kxky[ikxneg,ikyneg])

    # Now FT in y
    phi_kxy = np.zeros((nkx_no_negative, nky_inc_negative), dtype="complex")
    for ikx in range(0, nkx_no_negative):
        phi_kxy[ikx,:] = fft.ifft(phi_kxky_swap[ikx,:])

    # Now FT in x. Careful, because we want "c2r"; I think we can achieve this
    # by padding with the complex conjugate.
    phi_xy = np.zeros((nkx_inc_negative, nky_inc_negative), dtype="complex")
    for iky in range(0, nky_inc_negative):
        phi_kxy_padded = np.zeros((nkx_inc_negative), dtype="complex")
        phi_kxy_padded[:nkx_no_negative] = phi_kxy[:,iky]
        conjugate_vals = np.conj((phi_kxy[1:,iky])[::-1])
        #print("phi_kxy_padded, conjugate_vals = ", phi_kxy_padded, conjugate_vals)
        phi_kxy_padded[nkx_no_negative:] = conjugate_vals
        phi_xy[:,iky] = fft.ifft(phi_kxy_padded)


    # phi_kxky_full = np.zeros((nx, nky_full),dtype="complex")
    # # print("phi_kxky_full.shape = ", phi_kxky_full.shape)
    # phi_kxky_full[:,:nky] = phi_kxky
    # higher_entries = phi_kxky[:,1:]
    # # Reverse the order
    # print("higher_entries[0,:] = ", higher_entries[0,:])
    # higher_entries = higher_entries[:,::-1]
    # print("higher_entries[0,:] = ", higher_entries[0,:])
    # phi_kxky_full[:,nky:] = np.conj(higher_entries)

    #phi_xy = (fft.ifft2(phi_kxky_full))
    if np.max(abs(phi_xy.imag) > 1e-10):
        print("error! Complex part too big")
        print("phi_xy = ", phi_xy)
        sys.exit()
    # else:
    #     print("Success!")
    #     #sys.exit()
    return phi_xy.real

def make_phi2_t_movie(outnc_longname):
    """ """
    # Get phi(kx, ky, z, t)
    [t, kx, ky, z, phi2_kxky] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi2_vs_kxky")
    # print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
    # print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))

    z_idx = 1
    counter = 0
    # print("phi_t_ky_kx[0,:,:] = ", phi_t_ky_kx[0,:,:])
    # sys.exit()
    # For each time, perform the inverse Fourier transform to get phi(x,y). Make
    # this into a colorplot, and save it.
    print("kx = ", kx)
    print("ky = ", ky)
    dx = np.pi/np.max(kx)
    dy = np.pi/np.max(ky)
    nx = len(kx)
    ny = len(ky); nky=ny
    nky_full = nky*2-1  # length of ky including -ve freqs
    ny_full = nky_full
    xmax = (nx-1)*dx
    ymax = (ny_full-1)*dy
    xvals = np.linspace(0, xmax, nx)
    yvals = np.linspace(0, ymax, ny_full)
    # print("xvals = ", xvals)
    # print("yvals = ", yvals)
    xmesh, ymesh = np.meshgrid(xvals,yvals)

    # Want the colorbar to be the same for all time, and centered about zero.
    phi_xy_t = np.zeros((ny_full,nx,len(t)))

    for t_idx in range(0, len(t), 1):
        phi_kxky = phi2_kxky[t_idx, :,:]
        phi_xy = ifft_phi(kx, ky, phi_kxky)
        # print("phi_xy = ", phi_xy)
        # sys.exit()
        # Need to swap rows and columns; each column should be a single x
        phi_xy = phi_xy.T
        phi_xy_t[:,:,t_idx] = phi_xy

    max_phi = np.max(abs(phi_xy_t))
    print("max_phi = ", max_phi)

    for t_idx in range(0, len(t), 1):

        # x = fft.ifft(kx)
        # y = fft.ifft(ky)
        # print("x = ", x)
        # print("y = ", y)
        #print("phi_xy.shape = ", phi_xy.shape)
        #print("phi_xy = ", phi_xy)

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        # ax2 = fig.add_subplot(312, sharex=ax1)
        # ax3 = fig.add_subplot(313, sharex=ax1)
        pcm = ax1.pcolormesh(xvals, yvals, phi_xy_t[:,:,t_idx], cmap='inferno')#, vmin=-max_phi, vmax=max_phi)
        fig.colorbar(pcm, ax=ax1)
        #ax1.contourf(xmesh, ymesh, phi_xy)
        ax1.set_xlabel(r"$x$")
        ax1.set_ylabel(r"$y$")
        fig.suptitle("t={:.3f}".format(t[t_idx]))
        save_name="images/phi_t_{:02d}.png".format(counter)
        counter+=1
        #plt.show()
        # sys.exit()
        plt.savefig(save_name)
        plt.close()

    return

def make_phi_t_movie(outnc_longname, extra_upsample_fac=1):
    """ """
    # Get phi(kx, ky, z, t)
    [t, kx, ky, z, phi_vs_t] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi_vs_t")
    # print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
    # print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))

    z_idx = 1
    counter = 0
    phi_t_ky_kx = phi_vs_t[:, 0, z_idx, :, :]
    # print("phi_t_ky_kx[0,:,:] = ", phi_t_ky_kx[0,:,:])
    # sys.exit()
    # For each time, perform the inverse Fourier transform to get phi(x,y). Make
    # this into a colorplot, and save it.
    print("kx = ", kx)
    print("ky = ", ky)
    dx = np.pi/np.max(kx)
    dy = np.pi/np.max(ky)
    nx = int(((len(kx)-1)/2)*2*extra_upsample_fac+1)
    ny = len(ky); nky=ny
    nky_full = int((nky-1)*2*extra_upsample_fac+1)  # length of ky including -ve freqs
    ny_full = nky_full
    xmax = (nx-1)*dx
    ymax = (ny_full-1)*dy
    xvals = np.linspace(0, xmax, nx)
    yvals = np.linspace(0, ymax, ny_full)
    # print("xvals = ", xvals)
    # print("yvals = ", yvals)
    xmesh, ymesh = np.meshgrid(xvals,yvals)

    ## Want the colorbar to be the same for all time, and centered about zero.
    #phi_xy_t = np.zeros((ny_full,nx,len(t)))

    for t_idx in range(0, len(t), 1):
        phi_kxky = phi_t_ky_kx[t_idx, :,:]
        phi_xy = ifft_phi(kx, ky, phi_kxky, extra_upsample_fac=extra_upsample_fac)
        # print("phi_xy = ", phi_xy)
        # sys.exit()
        # Need to swap rows and columns; each column should be a single x
        phi_xy = phi_xy.T
        #phi_xy_t[:,:,t_idx] = phi_xy

    # max_phi = np.max(abs(phi_xy_t))
    # print("max_phi = ", max_phi)
    #
    # for t_idx in range(0, len(t), 1):

        # x = fft.ifft(kx)
        # y = fft.ifft(ky)
        # print("x = ", x)
        # print("y = ", y)
        #print("phi_xy.shape = ", phi_xy.shape)
        #print("phi_xy = ", phi_xy)

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        # ax2 = fig.add_subplot(312, sharex=ax1)
        # ax3 = fig.add_subplot(313, sharex=ax1)
        pcm = ax1.pcolormesh(xvals, yvals, abs(phi_xy), cmap='PRGn')#, vmin=-max_phi, vmax=max_phi)
        fig.colorbar(pcm, ax=ax1)
        #ax1.contourf(xmesh, ymesh, phi_xy)
        ax1.set_xlabel(r"$x$")
        ax1.set_ylabel(r"$y$")
        fig.suptitle("t={:.3f}".format(t[t_idx]))
        save_name="images/phi_t_{:02d}.png".format(counter)
        counter+=1
        #plt.show()
        # sys.exit()
        plt.savefig(save_name)
        plt.close()

    return

def make_phi_z_t_movie(outnc_longname, extra_upsample_fac=1, t_idxs=None):
    """ """
    # Get phi(kx, ky, z, t)
    [t, kx, ky, z, phi_vs_t] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi_vs_t")
    # print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
    # print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))

    counter = 0
    phi_t_z_ky_kx = phi_vs_t[:, 0, :, :, :]
    # print("phi_t_ky_kx[0,:,:] = ", phi_t_ky_kx[0,:,:])
    # sys.exit()
    # For each time, perform the inverse Fourier transform to get phi(x,y). Make
    # this into a colorplot, and save it.
    print("kx = ", kx)
    print("ky = ", ky)
    dx = np.pi/np.max(kx)
    dy = np.pi/np.max(ky)
    nx = int(((len(kx)-1)/2)*2*extra_upsample_fac+1)
    ny = len(ky); nky=ny
    nky_full = int((nky-1)*2*extra_upsample_fac+1)  # length of ky including -ve freqs
    ny_full = nky_full
    xmax = (nx-1)*dx
    ymax = (ny_full-1)*dy
    xvals = np.linspace(0, xmax, nx)
    yvals = np.linspace(0, ymax, ny_full)
    # print("xvals = ", xvals)
    # print("yvals = ", yvals)
    xmesh, ymesh = np.meshgrid(xvals,yvals)

    ## Want the colorbar to be the same for all time, and centered about zero.
    #phi_xy_t = np.zeros((ny_full,nx,len(t)))

    phi_xy_z = np.zeros((ny_full,nx,len(z)))
    if t_idxs is None:
        t_idxs = range(0, len(t), 20)
    for t_idx in t_idxs:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        phi_zkxky = phi_t_z_ky_kx[t_idx, :, :,:]
        for z_idx in range(0, len(z)):
            phi_xy = ifft_phi(kx, ky, phi_zkxky[z_idx,:,:])
            # Need to swap rows and columns; each column should be a single x
            phi_xy = phi_xy.T
            phi_xy_z[:,:,z_idx] = phi_xy

        for xidx in range(0, len(xvals)):
            for yidx in range(0, len(yvals)):
                ax1.plot(z, abs(phi_xy_z[yidx, xidx,:]), c="red", alpha=0.1)
        #ax1.contourf(xmesh, ymesh, phi_xy)
        ax1.set_xlabel(r"$z$")
        ax1.set_ylabel(r"$\vert \phi \vert$")
        fig.suptitle("t={:.3f}".format(t[t_idx]))
        save_name="images/phi_t_{:02d}.png".format(counter)
        counter+=1
        #plt.show()
        # sys.exit()
        plt.savefig(save_name)
        plt.close()

    return

def make_phi_z_movie(outnc_longname, extra_upsample_fac=1):
    """ """
    # Get phi(kx, ky, z, t)
    [t, kx, ky, z, phi_vs_t] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi_vs_t")
    # print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
    # print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))

    t_idx = -1
    counter = 0
    phi_z_ky_kx = phi_vs_t[t_idx, 0, :, :, :]
    # print("phi_t_ky_kx[0,:,:] = ", phi_t_ky_kx[0,:,:])
    # sys.exit()
    # For each time, perform the inverse Fourier transform to get phi(x,y). Make
    # this into a colorplot, and save it.
    print("kx = ", kx)
    print("ky = ", ky)
    dx = np.pi/np.max(kx)
    dy = np.pi/np.max(ky)
    nx = int(((len(kx)-1)/2)*2*extra_upsample_fac+1)
    ny = len(ky); nky=ny
    nky_full = int((nky-1)*2*extra_upsample_fac+1)  # length of ky including -ve freqs
    ny_full = nky_full
    xmax = (nx-1)*dx
    ymax = (ny_full-1)*dy
    xvals = np.linspace(0, xmax, nx)
    yvals = np.linspace(0, ymax, ny_full)
    # print("xvals = ", xvals)
    # print("yvals = ", yvals)
    xmesh, ymesh = np.meshgrid(xvals,yvals)

    ## Want the colorbar to be the same for all time, and centered about zero.
    #phi_xy_t = np.zeros((ny_full,nx,len(t)))

    for z_idx in range(0, len(z), 1):
        phi_kxky = phi_z_ky_kx[z_idx, :,:]
        phi_xy = ifft_phi(kx, ky, phi_kxky, extra_upsample_fac=extra_upsample_fac)
        # print("phi_xy = ", phi_xy)
        # sys.exit()
        # Need to swap rows and columns; each column should be a single x
        phi_xy = phi_xy.T
        #phi_xy_t[:,:,t_idx] = phi_xy

    # max_phi = np.max(abs(phi_xy_t))
    # print("max_phi = ", max_phi)
    #
    # for t_idx in range(0, len(t), 1):

        # x = fft.ifft(kx)
        # y = fft.ifft(ky)
        # print("x = ", x)
        # print("y = ", y)
        #print("phi_xy.shape = ", phi_xy.shape)
        #print("phi_xy = ", phi_xy)

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        # ax2 = fig.add_subplot(312, sharex=ax1)
        # ax3 = fig.add_subplot(313, sharex=ax1)
        pcm = ax1.pcolormesh(xvals, yvals, abs(phi_xy), cmap='inferno')#, vmin=-max_phi, vmax=max_phi)
        fig.colorbar(pcm, ax=ax1)
        #ax1.contourf(xmesh, ymesh, phi_xy)
        ax1.set_xlabel(r"$x$")
        ax1.set_ylabel(r"$y$")
        fig.suptitle("z={:.3f}".format(z[z_idx]))
        save_name="images/phi_z_{:02d}.png".format(counter)
        counter+=1
        #plt.show()
        # sys.exit()
        plt.savefig(save_name)
        plt.close()

    return

def make_phi_xyz_plots(outnc_longname):
    """ """
    # Get phi(kx, ky, z, t)
    [t, kx, ky, z, phi_vs_t] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi_vs_t")
    # print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
    # print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))

    t_idx = -1
    counter = 0
    phi_z_ky_kx = phi_vs_t[t_idx, 0, :, :, :]
    # print("phi_t_ky_kx[0,:,:] = ", phi_t_ky_kx[0,:,:])
    # sys.exit()
    # For each time, perform the inverse Fourier transform to get phi(x,y). Make
    # this into a colorplot, and save it.
    print("kx = ", kx)
    print("ky = ", ky)
    dx = np.pi/np.max(kx)
    dy = np.pi/np.max(ky)
    nx = len(kx)
    ny = len(ky); nky=ny
    nky_full = int((nky-1)*2+1)  # length of ky including -ve freqs
    ny_full = nky_full
    xmax = (nx-1)*dx
    ymax = (ny_full-1)*dy
    xvals = np.linspace(0, xmax, nx)
    yvals = np.linspace(0, ymax, ny_full)
    # print("xvals = ", xvals)
    # print("yvals = ", yvals)
    xmesh, ymesh = np.meshgrid(xvals,yvals)

    ## Want the colorbar to be the same for all time, and centered about zero.
    phi_xy_z = np.zeros((ny_full,nx,len(z)))

    for z_idx in range(0, len(z), 1):
        phi_kxky = phi_z_ky_kx[z_idx, :,:]
        phi_xy = ifft_phi(kx, ky, phi_kxky)
        # print("phi_xy = ", phi_xy)
        # sys.exit()
        # Need to swap rows and columns; each column should be a single x
        phi_xy = phi_xy.T
        phi_xy_z[:,:,z_idx] = phi_xy

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    for xidx in range(0, len(xvals)):
        for yidx in range(0, len(yvals)):
            ax1.plot(z, phi_xy_z[yidx, xidx,:], c="red", alpha=0.1)

    plt.show()

    return

def show_final_phi2_pic(outnc_longname):
    """ """
    # Get phi(kx, ky, z, t)
    # view_ncdf_variables(outnc_longname)
    [t, kx, ky, z, phi2] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi2_vs_kxky")
    # print("phi2.shape = ", phi2.shape)  # time, tube (?), z, kx, ky
    # print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))

    print("kx = ", kx)
    print("ky = ", ky)
    dx = np.pi/np.max(kx)
    dy = np.pi/np.max(ky)
    nx = len(kx)
    ny = len(ky); nky=ny
    nky_full = nky*2-1  # length of ky including -ve freqs
    ny_full = nky_full
    xmax = (nx-1)*dx
    ymax = (ny_full-1)*dy
    xvals = np.linspace(0, xmax, nx)
    yvals = np.linspace(0, ymax, ny_full)
    # print("xvals = ", xvals)
    # print("yvals = ", yvals)
    xmesh, ymesh = np.meshgrid(xvals,yvals)

    # Want the colorbar to be the same for all time, and centered about zero.

    phi2_kxky = phi2[-1, :,:]
    phi2_xy = ifft_phi(kx, ky, phi2_kxky)
    # print("phi_xy = ", phi_xy)
    # sys.exit()
    # Need to swap rows and columns; each column should be a single x
    phi2_xy = phi2_xy.T


    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # ax2 = fig.add_subplot(312, sharex=ax1)
    # ax3 = fig.add_subplot(313, sharex=ax1)
    pcm = ax1.pcolormesh(xvals, yvals, phi2_xy, cmap='RdBu_r') #vmin=-max_phi, vmax=max_phi)
    fig.colorbar(pcm, ax=ax1)
    #ax1.contourf(xmesh, ymesh, phi_xy)
    ax1.set_xlabel(r"$x$")
    ax1.set_ylabel(r"$y$")
    plt.show()

    return

def show_final_phi_pic(outnc_longname):
    """ """
    # Get phi(kx, ky, z, t)
    [t, kx, ky, z, phi_vs_t] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi_vs_t")
    print("phi2.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
    print("len(kx, ky, t, z) = ", len(kx), len(ky), len(t), len(z))

    z_idx = 1
    counter = 0
    phi_t_ky_kx = phi_vs_t[:, 0, z_idx, :, :]
    # print("phi_t_ky_kx[0,:,:] = ", phi_t_ky_kx[0,:,:])
    # sys.exit()
    # For each time, perform the inverse Fourier transform to get phi(x,y). Make
    # this into a colorplot, and save it.
    print("kx = ", kx)
    print("ky = ", ky)
    dx = np.pi/np.max(kx)
    dy = np.pi/np.max(ky)
    nx = len(kx)
    ny = len(ky); nky=ny
    nky_full = nky*2-1  # length of ky including -ve freqs
    ny_full = nky_full
    xmax = (nx-1)*dx
    ymax = (ny_full-1)*dy
    xvals = np.linspace(0, xmax, nx)
    yvals = np.linspace(0, ymax, ny_full)
    # print("xvals = ", xvals)
    # print("yvals = ", yvals)
    xmesh, ymesh = np.meshgrid(xvals,yvals)

    # Want the colorbar to be the same for all time, and centered about zero.

    phi_kxky = phi_t_ky_kx[-1, :,:]
    phi_xy = ifft_phi(kx, ky, phi_kxky)
    # print("phi_xy = ", phi_xy)
    # sys.exit()
    # Need to swap rows and columns; each column should be a single x
    phi_xy = phi_xy.T

    # x = fft.ifft(kx)
    # y = fft.ifft(ky)
    # print("x = ", x)
    # print("y = ", y)
    #print("phi_xy.shape = ", phi_xy.shape)
    #print("phi_xy = ", phi_xy)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # ax2 = fig.add_subplot(312, sharex=ax1)
    # ax3 = fig.add_subplot(313, sharex=ax1)
    pcm = ax1.pcolormesh(xvals, yvals, phi_xy, cmap='RdBu_r') #vmin=-max_phi, vmax=max_phi)
    fig.colorbar(pcm, ax=ax1)
    #ax1.contourf(xmesh, ymesh, phi_xy)
    ax1.set_xlabel(r"$x$")
    ax1.set_ylabel(r"$y$")
    plt.show()

    return

def compare_sims_with_different_nwrite():
    """ """
    nwrite100 = "example_nisl_nonlinear_only_vexb10_for_visualisation.out.nc"
    nwrite51 = "example_nisl_nonlinear_only_vexb10_for_visualisation_nwrite51.out.nc"
    nwrite25 = "example_nisl_nonlinear_only_vexb10_for_visualisation_nwrite25.out.nc"
    nwrite1 = "example_nisl_nonlinear_only_vexb10_for_visualisation_nwrite1.out.nc"
    nwrite1_smallerdt = "example_nisl_nonlinear_only_vexb10_for_visualisation_nwrite1_delt0.015.out.nc"

    [t_nwrite100, kx_nwrite100, ky_nwrite100, z_nwrite100, phi_vs_t_nwrite100] = extract_data_from_ncdf(nwrite100, "t", 'kx', 'ky', "zed", "phi_vs_t")
    [t_nwrite51, kx_nwrite51, ky_nwrite51, z_nwrite51, phi_vs_t_nwrite51] = extract_data_from_ncdf(nwrite51, "t", 'kx', 'ky', "zed", "phi_vs_t")
    [t_nwrite25, kx_nwrite25, ky_nwrite25, z_nwrite25, phi_vs_t_nwrite25] = extract_data_from_ncdf(nwrite25, "t", 'kx', 'ky', "zed", "phi_vs_t")
    [t_nwrite1, kx_nwrite1, ky_nwrite1, z_nwrite1, phi_vs_t_nwrite1] = extract_data_from_ncdf(nwrite1, "t", 'kx', 'ky', "zed", "phi_vs_t")
    [t_nwrite1_smallerdt, kx_nwrite1_smallerdt, ky_nwrite1_smallerdt,
        z_nwrite1_smallerdt, phi_vs_t_nwrite1_smallerdt] = extract_data_from_ncdf(nwrite1_smallerdt, "t", 'kx', 'ky', "zed", "phi_vs_t")

    # Plot |phi(t)| for each z, kx, ky
    nz_per_mode = len(z_nwrite51)

    z_idxs = [int(nz_per_mode/2)]  # The z idx of z=0
    if z_nwrite51[z_idxs[0]] != 0:
        print("Error! z[z_idx] = ", z[z_idx])
        sys.exit()
    # z_idx = 1 # Old
    counter = 0
    for ky_idx in range(0, len(ky_nwrite51)):
        for kx_idx in range(0, len(kx_nwrite51)):

            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            # ax2 = fig.add_subplot(312, sharex=ax1)
            # ax3 = fig.add_subplot(313, sharex=ax1)
            #for z_idx in range(0, nz_per_mode):
            ## For now, only plot for z=0

            for z_idx in z_idxs:
                # ax1.plot(t_nwrite100, abs(phi_vs_t_nwrite100[:, 0, z_idx, kx_idx,ky_idx]))
                # ax1.plot(t_nwrite51, abs(phi_vs_t_nwrite51[:, 0, z_idx, kx_idx,ky_idx]), ls="--")
                # ax1.plot(t_nwrite25, abs(phi_vs_t_nwrite25[:, 0, z_idx, kx_idx,ky_idx]), ls="-.")
                ax1.plot(t_nwrite1, abs(phi_vs_t_nwrite1[:, 0, z_idx, kx_idx,ky_idx]), ls="-", c="black", label="dt=0.03")
                ax1.plot(t_nwrite1_smallerdt, abs(phi_vs_t_nwrite1_smallerdt[:, 0, z_idx, kx_idx,ky_idx]), ls="-", c="red", label="dt=0.015")
                # ax1.plot(t_nwrite1, (phi_vs_t_nwrite1[:, 0, z_idx, kx_idx,ky_idx]).real, label="real(phi)")
                # ax1.plot(t_nwrite1, (phi_vs_t_nwrite1[:, 0, z_idx, kx_idx,ky_idx]).imag, label="im(phi)")
            ax1.set_xlabel(r"$t$")
            ax1.set_ylabel(r"$\phi$")
            ax1.legend(loc="best")
            fig.suptitle("kx={:.3f}, ky={:.3f}".format(kx_nwrite51[kx_idx], ky_nwrite51[ky_idx]))
            save_name="images/phi_t_{:02d}.png".format(counter)
            counter+=1
            #plt.savefig(save_name)
            #plt.close()
            plt.show()

    return

def plot_phi2t_for_folder(folder_name, nwrite):
    """For a folder of sims, plot phi2(t)
    Could also try to plot something like |G| (kxmax*U*delt) """
    filenames = glob.glob(folder_name + "/*.out.nc")

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    delt_vals = []
    amp_factor_vals = []
    amp_factor_vals2 = []
    amp_factor_vals3 = []
    amp_factor_vals4 = []
    for outnc_longname in filenames:
        [t, kx, ky, z, phi_vs_t] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi_vs_t")
        # print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, ky
        # print("kx = ", kx)  #   kx =  [ 0.         0.3334277  0.6668554  1.0002831 -1.0002831 -0.6668554
        #                     #           -0.3334277]
        # print("ky = ", ky)  #   ky =  [0.         0.33333333 0.66666667 1.         1.33333333]
        if len(t) > 1 :
            delt = (t[1] - t[0]) / nwrite
            #print("t = ", t)
            #print("delt = ", delt)
            # Look at kx=1, ky=1
            phi = phi_vs_t[:,0,0,4,4]
            if np.max(abs(phi)) < 1E4:
                ax1.plot(t, abs(phi), label="delt={:.6f}".format(delt))

            # Find the amplification factor |G| per timestep
            # First find the "good" bit of phi, where |phi| > 1E-15 but |phi| != NaN
            abs_phi = abs(phi)
            good_abs_phi = abs_phi[abs_phi<1E50]
            good_abs_phi = good_abs_phi[good_abs_phi>1E-10]
            if len(good_abs_phi) > 1:
                #print("delt, (good_abs_phi[-1]/good_abs_phi[0] = ", delt, (good_abs_phi[-1]/good_abs_phi[0]))
                delt_vals.append(delt)

                amp_factor_vals.append((good_abs_phi[-1]/good_abs_phi[-2])**(1/10))
                amp_factor_vals2.append((good_abs_phi[-2]/good_abs_phi[-3])**(1/10))
                amp_factor_vals3.append((good_abs_phi[-3]/good_abs_phi[-4])**(1/10))
                amp_factor_vals4.append((good_abs_phi[-1]/good_abs_phi[1])**(1/(10*len(good_abs_phi))))
            #print("min, max = ", np.min(good_abs_phi), np.max(good_abs_phi))

    #sys.exit()
    ax1.grid(True)
    ax1.legend(loc="best")
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_xlabel(r"$t$")
    ax1.set_ylabel(r"$\vert \phi \vert$")
    plt.show()


    print("delt_vals = ", delt_vals)
    print("amp_factor_vals = ", amp_factor_vals)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    marker_size = 50
    ax1.scatter(delt_vals, amp_factor_vals, marker="x", s=marker_size, label="sampling 10 points (1)")
    ax1.scatter(delt_vals, amp_factor_vals2, marker="x", s=marker_size, label="sampling 10 points (2)")
    ax1.scatter(delt_vals, amp_factor_vals3, marker="x", s=marker_size, label="sampling 10 points (3)")
    ax1.scatter(delt_vals, amp_factor_vals4, marker=".", s=200, label="sampling many points")
    #ax1.set_yscale("log")
    ax1.legend(loc="best")
    ax1.set_xlabel(r"$\beta \equiv k_x U_x \Delta t$")
    ax1.set_ylabel(r"$\vert G \vert$")
    ax1.grid(True)
    plt.show()

    return

def plot_phi2t_for_rk3_folder(folder_name, nwrite):
    """For a folder of sims, plot phi2(t)
    Could also try to plot something like |G| (kxmax*U*delt) """
    filenames = glob.glob(folder_name + "/*.out.nc")

    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    delt_vals = []
    amp_factor_vals = []
    for outnc_longname in filenames:
        [t, kx, ky, z, phi_vs_t] = extract_data_from_ncdf(outnc_longname, "t", 'kx', 'ky', "zed", "phi_vs_t")
        # print("phi_vs_t.shape = ", phi_vs_t.shape)  # time, tube (?), z, kx, k
        # print("kx = ", kx)  #   kx =  [ 0.         0.3334277  0.6668554  1.0002831 -1.0002831 -0.6668554
        #                     #           -0.3334277]
        # print("ky = ", ky)  #   ky =  [0.         0.33333333 0.66666667 1.         1.33333333]
        delt = (t[1] - t[0]) / nwrite
        # Look at kx=1, ky=1
        phi = phi_vs_t[:,0,0,4,4]
        # if np.max(abs(phi)) < 1E4:
        #     ax1.plot(t, abs(phi), label="delt={:.6f}".format(delt))

        # Find the amplification factor |G| per timestep
        # First find the "good" bit of phi, where |phi| > 1E-15 but |phi| != NaN
        abs_phi = abs(phi)
        good_abs_phi = abs_phi[abs_phi<1E50]
        good_abs_phi = good_abs_phi[good_abs_phi>1E-10]
        if len(good_abs_phi) > 1:
            #print("delt, (good_abs_phi[-1]/good_abs_phi[0] = ", delt, (good_abs_phi[-1]/good_abs_phi[0]))
            delt_vals.append(delt)

            amp_factor_vals.append((good_abs_phi[-1]/good_abs_phi[-2])**(1/10))
        #print("min, max = ", np.min(good_abs_phi), np.max(good_abs_phi))


    # ax1.grid(True)
    # ax1.legend(loc="best")
    # ax1.set_yscale("log")
    # ax1.set_xscale("log")
    # ax1.set_xlabel(r"$t$")
    # ax1.set_ylabel(r"$\vert \phi \vert$")
    # plt.show()

    print("delt_vals = ", delt_vals)
    print("amp_factor_vals = ", amp_factor_vals)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(delt_vals, amp_factor_vals, c="black", marker="x", label="RK3 sims")

    beta_min = np.min(delt_vals)
    beta_max = np.max(delt_vals)
    beta_vals = np.linspace(beta_min, beta_max, 1000)
    amp_factor_analytic = np.sqrt(1 - beta_vals**4 /12 + beta_vals**6/36)
    ax1.plot(beta_vals, amp_factor_analytic, label=r"analytic; $G = 1 - \beta^4/12 + \beta^6/36$")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$\beta \equiv k_x U_x \Delta t$")
    ax1.set_ylabel(r"$\vert G \vert$")
    ax1.legend(loc="best")
    ax1.grid(True)
    plt.show()

    return

if __name__ == "__main__":
    print("Hello world")
