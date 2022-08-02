"""Make plots for electromagnetic slab sims"""

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
# import ast

IMAGE_DIR = "./images/"

# basic_em_sim = "stella_sims/input_slab_ky0.1_explicit"
# mandell_beta1_kperp1 = "stella_sims/input_slab_ky0.1_explicit_mandell1"
# mandell_beta1_kperp001_new = "mandell_sims/input_slab_ky0.01_beta1_new"
# mandell_beta1_kperp1_new = "mandell_sims/input_slab_ky1_beta1_new"
# mandell_sf1_kperp1_new = "mandell_sims/input_slab_ky1_sf1_new"
# mandell_sf1_kperp001_new = "mandell_sims/input_slab_ky0.01_sf1_new"
# mandell_sf0002_kperp001_new = "mandell_sims/input_slab_ky0.01_sf0.002_new"
# mandell_sf00002_kperp001_new = "mandell_sims/input_slab_ky0.01_sf0.0002_new"
# mandell_sf000002_kperp001_new = "mandell_sims/input_slab_ky0.01_sf0.00002_new"
# mandell_sf10_kperp001_new = "mandell_sims/input_slab_ky0.01_sf10_new"
# mandell_beta1_kperp1_long_t = "stella_sims/input_slab_ky0.1_explicit_mandell2"
# mandell_beta1_kperp1_long_t_marconi = "mandell_sims/input_slab_ky0.1_beta1"

def examine_sim_output(sim_longname):
    """ """
    outnc_longname = sim_longname + ".out.nc"
    ### Plot geometric terms
    ## Code to compare geometry between stella and gs2

    # view_ncdf_variables(outnc_longname)
    # ['code_info', 'nproc', 'nmesh', 'ntubes', 'nkx', 'nky', 'nzed_tot',
    # 'nspecies', 'nmu', 'nvpa_tot', 't', 'charge', 'mass', 'dens', 'temp',
    # 'tprim', 'fprim', 'vnew', 'type_of_species', 'theta0', 'kx', 'ky', 'mu',
    # 'vpa', 'zed', 'bmag', 'gradpar', 'gbdrift', 'gbdrift0', 'cvdrift',
    # 'cvdrift0', 'kperp2', 'gds2', 'gds21', 'gds22', 'grho', 'jacob', 'q',
    # 'beta', 'shat', 'jtwist', 'drhodpsi', 'phi2', 'phi_vs_t', 'bpar2',
    # 'gvmus', 'gzvs', 'input_file']

    ###### Plot geometry
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    z, gds2, gds21, gds22, bmag, gradpar = extract_data_from_ncdf(outnc_longname,
                                    'zed', 'gds2', 'gds21', 'gds22', 'bmag', 'gradpar')
    ax1.plot(z, gds2)
    ax1.plot(z, gds21)
    ax1.scatter(z, gds22)
    ax2.plot(z, bmag)
    ax2.plot(z, gradpar)
    plt.show()

    t, z, phi_vs_t = extract_data_from_ncdf(outnc_longname, "t", 'zed', 'phi_vs_t')
    # print("len(t) = ", len(t))
    # print("len(z) = ", len(z))
    # print("phi_vs_t.shape = ", phi_vs_t.shape)  # [n_time, 1 , n_z, 1, 1]
    phi_vs_t = phi_vs_t[:,0,:,0,0]
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    ax1.plot(z, phi_vs_t[-20,:].real, label="t[-20]")
    ax2.plot(z, phi_vs_t[-20,:].imag)
    ax1.plot(z, phi_vs_t[-10,:].real, label="t[-10]")
    ax2.plot(z, phi_vs_t[-10,:].imag)
    ax1.plot(z, phi_vs_t[-1,:].real, label="t[-1]")
    ax2.plot(z, phi_vs_t[-1,:].imag)
    ax2.set_xlabel("z")
    ax2.set_ylabel("Im(phi)")
    ax1.set_ylabel("Re(phi)")
    ax1.legend(loc="best")
    for ax in [ax1, ax2]:
        ax.grid(True)
    plt.show()

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    ax1.plot(t, phi_vs_t[:,int(len(z)*0.5)].real, label=("z=" + str(z[int(len(z)*0.5)]) ))
    ax2.plot(t, phi_vs_t[:,int(len(z)*0.5)].imag)
    ax1.plot(t, phi_vs_t[:,int(len(z)*0.2)].real, label=("z=-1.9" ))
    ax2.plot(t, phi_vs_t[:,int(len(z)*0.2)].imag)
    ax2.set_xlabel("t")
    ax2.set_ylabel("Im(phi)")
    ax1.set_ylabel("Re(phi)")
    ax1.legend(loc="best")
    for ax in [ax1, ax2]:
        ax.grid(True)
    plt.show()

    t_chop = int(len(t)/2)
    print("len(t), t_chop, t[t_chop] = ", len(t), t_chop, t[t_chop] )
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    ax1.plot(t[:t_chop], phi_vs_t[:t_chop,int(len(z)*0.5)].real, label=("z=" + str(z[int(len(z)*0.5)]) ))
    ax2.plot(t[:t_chop], phi_vs_t[:t_chop,int(len(z)*0.5)].imag)
    ax1.plot(t[:t_chop], phi_vs_t[:t_chop,int(len(z)*0.2)].real, label=("z=-1.9" ))
    ax2.plot(t[:t_chop], phi_vs_t[:t_chop,int(len(z)*0.2)].imag)
    ax2.set_xlabel("t")
    ax2.set_ylabel("Im(phi)")
    ax1.set_ylabel("Re(phi)")
    ax1.legend(loc="best")
    for ax in [ax1, ax2]:
        ax.grid(True)
    plt.show()

    return

def examine_first_sim():
    """ """
    examine_sim_output(basic_em_sim)
    return

def find_ksaw_properties_from_pickle(phi_vs_t_file):
    """ """
    file = open(phi_vs_t_file, "rb")
    [z, t, phiz_final, phit_mid] = pickle.load(file)
    file.close()
    ### Lines are :
    #  "z"
    #  z
    #  "t"
    #  t
    #  "phi_vs_t[-1,:]"
    #  phi_vs_t[-1,:]
    #  "phi_vs_t[:,int(len(z)*0.5)]"
    #  phi_vs_t[:,int(len(z)*0.5)]

    # z = np.array(ast.literal_eval(lines[1]))
    # t = np.array(ast.literal_eval(lines[3]))
    # phiz_final = ast.literal_eval(lines[5])
    # print("phiz_final = ", phiz_final)
    # sys.exit()
    # phit_mid = np.array(ast.literal_eval(lines[7]))

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    ax1.plot(z, phiz_final.real)
    ax2.plot(z, phiz_final.imag)
    ax2.set_xlabel("z")
    ax2.set_ylabel("Im(phi)")
    ax1.set_ylabel("Re(phi)")
    fig.suptitle(phi_vs_t_file)
    for ax in [ax1, ax2]:
        ax.grid(True)
    plt.show()

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    ax1.plot(t, phit_mid.real)
    ax2.plot(t, phit_mid.imag)
    ax2.set_xlabel("t")
    ax2.set_ylabel("Im(phi)")
    ax1.set_ylabel("Re(phi)")
    ax1.legend(loc="best")
    fig.suptitle(phi_vs_t_file)
    for ax in [ax1, ax2]:
        ax.grid(True)
    plt.show()


    return

def find_ksaw_properties_from_outnc(outnc_longname):
    """ """
    t, z, phi_vs_t, beta = extract_data_from_ncdf(outnc_longname, "t", 'zed', 'phi_vs_t', 'beta')
    print("beta=",beta)
    # print("len(t) = ", len(t))
    # print("len(z) = ", len(z))
    # print("phi_vs_t.shape = ", phi_vs_t.shape)  # [n_time, 1 , n_z, 1, 1]
    phi_vs_t = phi_vs_t[:,0,:,0,0]
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    ax1.plot(z, phi_vs_t[0,:].real, label="t[0]")
    ax2.plot(z, phi_vs_t[0,:].imag)
    ax1.plot(z, phi_vs_t[1,:].real, label="t[1]")
    ax2.plot(z, phi_vs_t[1,:].imag)
    ax1.plot(z, phi_vs_t[-10,:].real, label="t[-10]")
    ax2.plot(z, phi_vs_t[-10,:].imag)
    ax1.plot(z, phi_vs_t[-1,:].real, label="t[-1]")
    ax2.plot(z, phi_vs_t[-1,:].imag)
    ax2.set_xlabel("z")
    ax2.set_ylabel("Im(phi)")
    ax1.set_ylabel("Re(phi)")
    ax1.legend(loc="best")
    for ax in [ax1, ax2]:
        ax.grid(True)
    plt.show()

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    ax1.plot(t, phi_vs_t[:,int(len(z)*0.5)].real, label=("z=" + str(z[int(len(z)*0.5)]) ))
    ax2.plot(t, phi_vs_t[:,int(len(z)*0.5)].imag)
    ax1.plot(t, phi_vs_t[:,int(len(z)*0.2)].real, label=("z=-1.9" ))
    ax2.plot(t, phi_vs_t[:,int(len(z)*0.2)].imag)
    ax2.set_xlabel("t")
    ax2.set_ylabel("Im(phi)")
    ax1.set_ylabel("Re(phi)")
    ax1.legend(loc="best")
    for ax in [ax1, ax2]:
        ax.grid(True)
    plt.show()

    return

def find_ksaw_properties_from_outnc_with_apar(outnc_longname):
    """ """
    t, z, phi_vs_t, apar_vs_t, beta = extract_data_from_ncdf(outnc_longname, "t", 'zed', 'phi_vs_t', 'apar_vs_t', 'beta')
    print("beta=",beta)
    # print("len(t) = ", len(t))
    # print("len(z) = ", len(z))
    # print("phi_vs_t.shape = ", phi_vs_t.shape)  # [n_time, 1 , n_z, 1, 1]
    phi_vs_t = phi_vs_t[:,0,:,0,0]
    apar_vs_t = apar_vs_t[:,0,:,0,0]
    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(223, sharex=ax1)
    ax3 = fig.add_subplot(222)
    ax4 = fig.add_subplot(224, sharex=ax3)
    ax1.plot(z, phi_vs_t[0,:].real, label="t[0]")
    ax2.plot(z, phi_vs_t[0,:].imag)
    ax1.plot(z, phi_vs_t[1,:].real, label="t[1]")
    ax2.plot(z, phi_vs_t[1,:].imag)
    ax1.plot(z, phi_vs_t[-10,:].real, label="t[-10]")
    ax2.plot(z, phi_vs_t[-10,:].imag)
    ax1.plot(z, phi_vs_t[-1,:].real, label="t[-1]")
    ax2.plot(z, phi_vs_t[-1,:].imag)
    ax3.plot(z, apar_vs_t[0,:].real, label="t[0]")
    ax4.plot(z, apar_vs_t[0,:].imag)
    ax3.plot(z, apar_vs_t[1,:].real, label="t[1]")
    ax4.plot(z, apar_vs_t[1,:].imag)
    ax3.plot(z, apar_vs_t[-10,:].real, label="t[-10]")
    ax4.plot(z, apar_vs_t[-10,:].imag)
    ax3.plot(z, apar_vs_t[-1,:].real, label="t[-1]")
    ax4.plot(z, apar_vs_t[-1,:].imag)
    ax2.set_xlabel("z")
    ax4.set_xlabel("z")
    ax2.set_ylabel("Im(phi)")
    ax1.set_ylabel("Re(phi)")
    ax4.set_ylabel("Im(apar)")
    ax3.set_ylabel("Re(apar)")
    ax1.legend(loc="best")
    for ax in [ax1, ax2, ax3, ax4]:
        ax.grid(True)
    plt.show()

    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(223, sharex=ax1)
    ax3 = fig.add_subplot(222)
    ax4 = fig.add_subplot(224, sharex=ax3)
    ax1.plot(t, phi_vs_t[:,int(len(z)*0.5)].real, label=("z=" + str(z[int(len(z)*0.5)]) ))
    ax2.plot(t, phi_vs_t[:,int(len(z)*0.5)].imag)
    ax1.plot(t, phi_vs_t[:,int(len(z)*0.2)].real, label=("z=-1.9" ))
    ax2.plot(t, phi_vs_t[:,int(len(z)*0.2)].imag)
    ax3.plot(t, apar_vs_t[:,int(len(z)*0.5)].real, label=("z=" + str(z[int(len(z)*0.5)]) ))
    ax4.plot(t, apar_vs_t[:,int(len(z)*0.5)].imag)
    ax3.plot(t, apar_vs_t[:,int(len(z)*0.2)].real, label=("z=-1.9" ))
    ax4.plot(t, apar_vs_t[:,int(len(z)*0.2)].imag)
    ax2.set_xlabel("t")
    ax4.set_xlabel("t")
    ax2.set_ylabel("Im(phi)")
    ax1.set_ylabel("Re(phi)")
    ax2.set_ylabel("Im(apar)")
    ax1.set_ylabel("Re(apar)")
    ax1.legend(loc="best")
    for ax in [ax1, ax2, ax3, ax4]:
        ax.grid(True)
    plt.show()

    return

def compare_sims(outnc_longnames, sim_types, labels, normalise=False, scatter=False,
                 title=None, include_im = False):
    """ """

    t_list = []
    z_list = []
    phi_vs_t_list = []
    linestyles=((0,(1,0)), (0, (2,2)), (0,(5,1,1,1)), (0,(1,1)), (0,(1,0)), (0, (2,2)))
    linewidths = [4, 4, 3, 3, 2, 2]
    for idx, outnc_longname in enumerate(outnc_longnames):
        if sim_types[idx] == "stella":
            t, z, phi_vs_t, beta = extract_data_from_ncdf_with_xarray(outnc_longname, "t", 'zed', 'phi_vs_t', 'beta')
            # print("stella. phi_vs_t.shape = ", phi_vs_t.shape)
            phi_vs_t = np.array(phi_vs_t[:,0,:,0,0,:])
            # shape is now (time, zed, ri)
            # print("stella. phi_vs_t.shape = ", phi_vs_t.shape)
        elif sim_types[idx] == "gs2":
            t, z, phi_vs_t, beta = extract_data_from_ncdf_with_xarray(outnc_longname, "t", 'theta', 'phi_t', 'beta')
            # print("gs2. phi_vs_t.shape = ", phi_vs_t.shape)
            phi_vs_t = np.array(phi_vs_t[:,0,0,:,:])
            # shape is now (time, zed, ri)
            # print("gs2. phi_vs_t.shape = ", phi_vs_t.shape)
        else:
            print("sim type not recognised! Aborting")
            sys.exit()
        if normalise:
            phi_vs_t = phi_vs_t/np.max(abs(phi_vs_t))
        t_list.append(t); z_list.append(z); phi_vs_t_list.append(phi_vs_t);
        # print("len(t) = ", len(t))
        # print("len(z) = ", len(z))
        # sys.exit()
    fig = plt.figure(figsize=(12,8))
    if include_im:
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex=ax1)
        lower_ax = ax2
        axes = [ax1, ax2]
    else:
        ax1 = fig.add_subplot(111)
        lower_ax = ax1
        axes = [ax1]
    for idx, outnc_longname in enumerate(outnc_longnames):
        # print("outnc_longname=", outnc_longname)
        t = t_list[idx]; z = z_list[idx]; phi_vs_t = phi_vs_t_list[idx];
        # print("t=", t)
        # print("phi_vs_t[:,int(len(z)*0.5),0]=", phi_vs_t[:,int(len(z)*0.5),0])
        ax1.plot(t, phi_vs_t[:,int(len(z)*0.5),0], label=(labels[idx] + ", z=" + str(float(z[int(len(z)*0.5)]))), ls=linestyles[idx], lw=linewidths[idx])
        if scatter:
            ax1.scatter(t, phi_vs_t[:,int(len(z)*0.5),0], marker="x", s=40)
        if include_im:
            ax2.plot(t, phi_vs_t[:,int(len(z)*0.5),1], ls=linestyles[idx], lw=linewidths[idx])
    lower_ax.set_xlabel("t")
    if include_im:
        ax2.set_ylabel("Im(phi)")
    ax1.set_ylabel("Re(phi)")
    ax1.legend(loc="best")
    for ax in axes:
        ax.grid(True)
    if title is not None:
        fig.suptitle(title)
    plt.tight_layout()
    plt.show()

    fig = plt.figure(figsize=(12,8))
    if include_im:
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex=ax1)
        lower_ax = ax2
        axes = [ax1, ax2]
    else:
        ax1 = fig.add_subplot(111)
        lower_ax = ax1
        axes = [ax1]
    for idx, outnc_longname in enumerate(outnc_longnames):
        t = t_list[idx]; z = z_list[idx]; phi_vs_t = phi_vs_t_list[idx];
        ax1.plot(z, phi_vs_t[0,:,0], label=(labels[idx] + ", t=0"), ls=linestyles[idx], lw=linewidths[idx])
        if include_im:
            ax2.plot(z, phi_vs_t[0,:,1], ls=linestyles[idx], lw=linewidths[idx])
    lower_ax.set_xlabel("z")
    if include_im:
        ax2.set_ylabel("Im(phi)")
    ax1.set_ylabel("Re(phi)")
    ax1.legend(loc="best")
    for ax in axes:
        ax.grid(True)
    if title is not None:
        fig.suptitle(title)
    plt.tight_layout()
    plt.show()

    fig = plt.figure(figsize=(12,8))
    if include_im:
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex=ax1)
        lower_ax = ax2
        axes = [ax1, ax2]
    else:
        ax1 = fig.add_subplot(111)
        lower_ax = ax1
        axes = [ax1]
    for idx, outnc_longname in enumerate(outnc_longnames):
        t = t_list[idx]; z = z_list[idx]; phi_vs_t = phi_vs_t_list[idx];
        ax1.plot(z, phi_vs_t[-1,:,0], label=(labels[idx] + ", t=tfinal"), ls=linestyles[idx], lw=linewidths[idx])
        if include_im:
            ax2.plot(z, phi_vs_t[-1,:,1], ls=linestyles[idx], lw=linewidths[idx])
    lower_ax.set_xlabel("z")
    if include_im:
        ax2.set_ylabel("Im(phi)")
    ax1.set_ylabel("Re(phi)")
    ax1.legend(loc="best")
    for ax in axes:
        ax.grid(True)
    if title is not None:
        fig.suptitle(title)
    plt.tight_layout()
    plt.show()

    return

def make_plots_for_movies(outnc_longname, sim_type, output_folder):
    """ """
    if sim_type == "stella":
        t, z, phi_vs_t, beta = extract_data_from_ncdf_with_xarray(outnc_longname, "t", 'zed', 'phi_vs_t', 'beta')
        phi_vs_t = np.array(phi_vs_t[:,0,:,0,0,:])
    elif sim_type == "gs2":
        t, z, phi_vs_t, beta = extract_data_from_ncdf_with_xarray(outnc_longname, "t", 'theta', 'phi_t', 'beta')
        # print("gs2. phi_vs_t.shape = ", phi_vs_t.shape)
        phi_vs_t = np.array(phi_vs_t[:,0,0,:,:])
    else:
        print("sim type not supporting, aborting")
        sys.exit()

    # print("stella. phi_vs_t.shape = ", phi_vs_t.shape)

    # Work out the max. vals, so they're the same every time
    min_z = np.max(z/np.pi)*1.01
    max_z = np.min(z/np.pi)*1.01
    min_t = np.min(t)*1.01
    max_t = np.max(t)*1.01
    max_real = np.max(phi_vs_t[:,:,0])*1.05
    min_real = np.min(phi_vs_t[:,:,0])*1.05
    max_imag = np.max(phi_vs_t[:,:,1])*1.05
    min_imag = np.min(phi_vs_t[:,:,1])*1.05
    print("max_real, min_real = ", max_real, min_real)
    if (max_imag - min_imag) < 1E-12:
        max_imag += 1E-12
        min_imag -= 1E-12
    if (max_real - min_real) < 1E-12:
        max_real += 1E-12
        min_real -= 1E-12
    print("max_real, min_real = ", max_real, min_real)

    for tidx, tval in enumerate(t):
        fig = plt.figure(figsize=(12,12))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex=ax1)
        ax1.plot(z/np.pi, phi_vs_t[tidx,:,0], label=(f't={float(tval):.3f}'), lw=3)
        ax2.plot(z/np.pi, phi_vs_t[tidx,:,1], lw=3)
        ax2.set_xlabel(r"$z/\pi$")
        ax2.set_ylabel("Im(phi)")
        ax1.set_ylabel("Re(phi)")
        ax1.legend(loc="best")
        for ax in [ax1, ax2]:
            ax.grid(True)
            ax.set_xlim([min_z, max_z])
        ax1.set_ylim([min_real, max_real])
        ax2.set_ylim([min_imag, max_imag])
        save_name = output_folder + f'/{tidx:03}'
        plt.tight_layout()
        plt.savefig(save_name)
        plt.close()

    fig = plt.figure(figsize=(12,12))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    for zidx in [0, int(len(z)/8), int(len(z)/4), int(3*len(z)/8), int(len(z)/2)]:
        zval = z[zidx]
        ax1.plot(t, phi_vs_t[:, zidx,0], lw=3)
        ax2.plot(t, phi_vs_t[:, zidx,1], label=(f'z/pi={float(zval/np.pi):.3f}'), lw=3)
    for zidx in [int(5*len(z)/8), int(3/4*len(z)), int(7*len(z)/8), -1]:
        zval = z[zidx]
        ax1.plot(t, phi_vs_t[:, zidx,0], lw=3, ls="--")
        ax2.plot(t, phi_vs_t[:, zidx,1], label=(f'z/pi={float(zval/np.pi):.3f}'), lw=3, ls="--")
    ax2.set_xlabel(r"$t$")
    ax2.set_ylabel("Im(phi)")
    ax1.set_ylabel("Re(phi)")
    ax2.legend(loc="best")
    for ax in [ax1, ax2]:
        ax.grid(True)
        ax.set_xlim([min_t, max_t])
    ax1.set_ylim([min_real, max_real])
    ax2.set_ylim([min_imag, max_imag])
    save_name = output_folder + '/phi_t.png'
    plt.tight_layout()
    plt.savefig(save_name)
    plt.close()
    return

def compare_sims_with_apar(outnc_longnames, labels):
    """ """

    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(223, sharex=ax1)
    ax3 = fig.add_subplot(222)
    ax4 = fig.add_subplot(224, sharex=ax3)
    for idx, outnc_longname in enumerate(outnc_longnames):
        t, z, phi_vs_t, apar_vs_t, beta = extract_data_from_ncdf(outnc_longname, "t", 'zed', 'phi_vs_t', 'apar_vs_t', 'beta')
        phi_vs_t = phi_vs_t[:,0,:,0,0]
        apar_vs_t = apar_vs_t[:,0,:,0,0]
        print("len(t) = ", len(t))
        print("len(z) = ", len(z))
        print("phi_vs_t.shape = ", phi_vs_t.shape)
        # print("phi_vs_t = ", phi_vs_t)
        # print("phi_vs_t.imag = ", phi_vs_t.imag)
        # sys.exit()
        ax1.plot(z, phi_vs_t[0,:].real, label=(labels[idx] + ", t=0"))
        ax2.plot(z, phi_vs_t[0,:].imag)
        ax3.plot(z, apar_vs_t[0,:].real, label=(labels[idx] + ", t=0" ))
        ax4.plot(z, apar_vs_t[0,:].imag)
    ax2.set_xlabel("z")
    ax4.set_xlabel("z")
    ax2.set_ylabel("Im(phi)")
    ax1.set_ylabel("Re(phi)")
    ax4.set_ylabel("Im(apar)")
    ax3.set_ylabel("Re(apar)")
    ax1.legend(loc="best")
    ax3.legend(loc="best")
    for ax in [ax1, ax2, ax3, ax4]:
        ax.grid(True)
    plt.show()

    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(223, sharex=ax1)
    ax3 = fig.add_subplot(222)
    ax4 = fig.add_subplot(224, sharex=ax3)
    for idx, outnc_longname in enumerate(outnc_longnames):
        t, z, phi_vs_t, apar_vs_t, beta = extract_data_from_ncdf(outnc_longname, "t", 'zed', 'phi_vs_t', 'apar_vs_t', 'beta')
        phi_vs_t = phi_vs_t[:,0,:,0,0]
        apar_vs_t = apar_vs_t[:,0,:,0,0]
        print("len(t) = ", len(t))
        print("len(z) = ", len(z))
        print("phi_vs_t.shape = ", phi_vs_t.shape)
        # print("phi_vs_t = ", phi_vs_t)
        # print("phi_vs_t.imag = ", phi_vs_t.imag)
        # sys.exit()
        ax1.plot(t, phi_vs_t[:,int(len(z)*0.5)].real, label=(labels[idx] + ", z=" + str(z[int(len(z)*0.5)]) ))
        ax2.plot(t, phi_vs_t[:,int(len(z)*0.5)].imag)
        ax3.plot(t, apar_vs_t[:,int(len(z)*0.25)].real, label=(labels[idx] + ", z=" + str(z[int(len(z)*0.25)]) ))
        ax4.plot(t, apar_vs_t[:,int(len(z)*0.25)].imag)
    ax2.set_xlabel("t")
    ax4.set_xlabel("t")
    ax2.set_ylabel("Im(phi)")
    ax1.set_ylabel("Re(phi)")
    ax4.set_ylabel("Im(apar)")
    ax3.set_ylabel("Re(apar)")
    ax1.legend(loc="best")
    ax3.legend(loc="best")
    for ax in [ax1, ax2, ax3, ax4]:
        ax.grid(True)
    plt.show()

    return

def benchmark_stella_vs_gs2():
    """Compare similar simulations between stella and GS2"""

    stella_sim = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar1/input_fully_explicit.out.nc"
    stella_sim_short = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar1/input_fully_explicit_short.out.nc"
    stella_sim_new_init = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar1/input_fully_explicit_init_matching_gs2.out.nc"
    stella_sim_new_init_ntheta144 = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar1/input_fully_explicit_init_matching_gs2_ntheta144.out.nc"
    gs2_long_sim = "sims/gs2_ky1_beta1_zero_gradients_fapar1_fbpar1/input_long.out.nc"
    compare_sims([gs2_long_sim, stella_sim, stella_sim_new_init, stella_sim_new_init_ntheta144],
            ["gs2", "stella", "stella", "stella"],
            ["gs2", "stella (explicit)",
             "stella (explicit, gs2-like init)", "stella (explicit, gs2-like init) nzed=144"], normalise=False,
             title="fapar=1, fbpar=1")

    # stella_sim = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar1/input_fully_explicit.out.nc"
    # stella_sim_ntheta72 = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar1/input_fully_explicit_ntheta72.out.nc"
    # stella_sim_ntheta144 = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar1/input_fully_explicit_ntheta144.out.nc"
    # gs2_long_sim = "sims/gs2_ky1_beta1_zero_gradients_fapar1_fbpar1/input_long.out.nc"
    # compare_sims([gs2_long_sim, stella_sim, stella_sim_ntheta72, stella_sim_ntheta144],
    #         ["gs2", "stella", "stella", "stella"],
    #         ["gs2", "stella (explicit) nzed=36", "stella (explicit) nzed=72", "stella (explicit) nzed=144"],
    #         normalise=False)
    # stella_sim = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar1/input_fully_explicit.out.nc"
    # stella_sim_new = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar1/input_fully_explicit_new.out.nc"
    # gs2_long_sim = "sims/gs2_ky1_beta1_zero_gradients_fapar1_fbpar1/input_long.out.nc"
    # compare_sims([gs2_long_sim, stella_sim, stella_sim_new],
    #         ["gs2", "stella", "stella"], ["gs2", "stella (explicit)", "stella (explicit) new"], normalise=False)
    #
    # stella_sim = "sims/stella_ky1_beta1_zero_gradients_fapar0_fbpar1/input_fully_explicit.out.nc"
    # gs2_long_sim = "sims/gs2_ky1_beta1_zero_gradients_fapar0_fbpar1/input_long.out.nc"
    # compare_sims([gs2_long_sim, stella_sim], ["gs2", "stella"], ["gs2", "stella (explicit)"], normalise=True)
    #


    # gs2_long_sim = "sims/gs2_ky1_beta1_zero_gradients_fapar0_fbpar0/input_long.out.nc"
    # stella_sim = "sims/stella_ky1_beta1_zero_gradients_fapar0_fbpar0/input_fully_explicit.out.nc"
    # stella_sim_short = "sims/stella_ky1_beta1_zero_gradients_fapar0_fbpar0/input_fully_explicit_short.out.nc"
    # stella_sim_str_impl = "sims/stella_ky1_beta1_zero_gradients_fapar0_fbpar0/input_implicit_streaming.out.nc"
    # compare_sims([gs2_long_sim, stella_sim, stella_sim_short, stella_sim_str_impl],
    #              ["gs2", "stella", "stella", "stella"],
    #              ["gs2", "stella (explicit)", "stella (explicit, short)", "stella (implicit)"],
    #              normalise=True)
    # stella_sim_nwrite100 = "sims/stella_ky1_beta1_zero_gradients_fapar0_fbpar0/input_fully_explicit_nwrite100.out.nc"
    # stella_sim_nwrite10 = "sims/stella_ky1_beta1_zero_gradients_fapar0_fbpar0/input_fully_explicit_short.out.nc"
    # stella_sim_nwrite25 = "sims/stella_ky1_beta1_zero_gradients_fapar0_fbpar0/input_fully_explicit_mid.out.nc"
    # compare_sims([stella_sim_nwrite10, stella_sim_nwrite25, stella_sim_nwrite100],
    #              ["stella", "stella", "stella"],
    #              ["nwrite=10", "nwrite=25", "nwrite=100"],
    #               normalise=False, scatter=True)
    return

def benchmark_stella_vs_gs2_fapar1_fbpar0():
    """ """
    gs2_sim = "sims/gs2_ky1_beta1_zero_gradients_fapar1_fbpar0/input_long.out.nc"
    gs2_sim_tup0 = "sims/gs2_ky1_beta1_zero_gradients_fapar1_fbpar0/input_long_fexpr_0.5.out.nc"
    gs2_sim_tup0_zup0 = "sims/gs2_ky1_beta1_zero_gradients_fapar1_fbpar0/input_long_fexpr_0.5_bakdif0.out.nc"
    stella_sim = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_fully_explicit.out.nc"
    stella_sim_str_impl = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_implicit_streaming.out.nc"
    stella_sim_str_impl_zup0_tup0 = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_implicit_streaming_zup0_tup0.out.nc"
    stella_sim_str_impl_zup0 = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_implicit_streaming_zup0.out.nc"
    stella_sim_str_impl_tup0 = "sims/stella_ky1_beta1_zero_gradients_fapar1_fbpar0/input_implicit_streaming_tup0.out.nc"
    compare_sims([gs2_sim, stella_sim, stella_sim_str_impl,
                  stella_sim_str_impl_zup0, stella_sim_str_impl_tup0,
                  stella_sim_str_impl_zup0_tup0],
                 ["gs2", "stella", "stella", "stella", "stella",
                 "stella", "stella",],
                 ["gs2", "stella (explicit)", "stella (implicit, zupw=tupw=0.02)",
                  "stella (implicit, zupw=0, tupw=0.02)", "stella (implicit, zupw=0.2, tupw=0)",
                  "stella (implicit, zupw=tupw=0)"],
                 normalise=False,
                 title="fapar=1, fbpar=0")

    compare_sims([gs2_sim, gs2_sim_tup0, gs2_sim_tup0_zup0,
                  stella_sim_str_impl,
                  stella_sim_str_impl_zup0_tup0],
                 ["gs2", "gs2", "gs2", "stella", "stella",
                  ],
                 ["gs2", "gs2 (tupw=0)", "gs2 (zupw=tupw=0)",
                  "stella (implicit, zupw=0, tupw=0.02)",
                  "stella (implicit, zupw=tupw=0)"],
                 normalise=False,
                 title="fapar=1, fbpar=0")

    return

def examine_second_sim():
    """ """
    examine_sim_output(mandell_beta1_kperp1_long_t_marconi)
    return

def examine_gs2_sim():
    """Take a look at an GS2 simulation, which has the following parameters:
        beta = 1
        aky=1
        ntheta=36
        nperiod=1
        boundary_option = "periodic"
        fapar=1
        fbpar=0
        nspec=2
        fprim_i=fprim_e = 0.809
        tprim_i=tprim_e=2.537
        ginit_option = "single_parallel_mode"
        ikpar_init = 1
        kpar_init = 1
        """
    print("Hello world")
    outnc_longname = "gs2_slab_sims/input.out.nc"
    view_ncdf_variables(outnc_longname)
    [shear, theta, t, phi, apar] = extract_data_from_ncdf(outnc_longname, "shat", "theta", "t", "phi_t", "apar_t")
    print("shear = ", shear)
    print("theta.shape = ", theta.shape)
    print("t.shape = ", t.shape)
    print("phi.shape = ", phi.shape)
    print("apar.shape = ", apar.shape)
    phi_for_theta_0 = phi[:,0,0,int(len(theta)/2)]
    phi_for_t_0 = phi[0,0,0,:]
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.plot(theta/np.pi, phi_for_t_0.real)
    ax1.plot(theta/np.pi, phi_for_t_0.imag)
    ax1.plot(theta/np.pi, abs(phi_for_t_0))
    ax2.plot(t, phi_for_theta_0.real)
    ax2.plot(t, phi_for_theta_0.imag)
    ax2.plot(t, abs(phi_for_theta_0))

    plt.show()
    return

if __name__ == "__main__":
    print("Hello world")
    # examine_first_sim()
    # examine_second_sim()
    # benchmark_stella_vs_mandell()
    # benchmark_stella_vs_mandell2()
    # find_ksaw_properties_from_outnc(mandell_beta1_kperp001_new + ".out.nc")
    # find_ksaw_properties_from_outnc(mandell_sf1_kperp001_new + ".out.nc")
    # find_ksaw_properties_from_outnc("mandell_sims/input_thursday1.out.nc")
    # find_ksaw_properties_from_outnc_with_apar("mandell_sims/input_thursday17_nions_beta1_write_apar.out.nc")
    # find_ksaw_properties_from_outnc_with_apar("mandell_sims/input_thursday18_wapar_upar.out.nc")
    # find_ksaw_properties_from_outnc_with_apar("mandell_sims/input_thursday19_wapar_upar_beta10.out.nc")

    # compare_sims([
    # #               #"mandell_sims/input_thursday2_dt.out.nc",
    # #               #"mandell_sims/input_thursday5_ky.out.nc",
    # #               "mandell_sims/input_thursday8_nions.out.nc",
    # #               "mandell_sims/input_thursday12_smions_dt.out.nc",
    #                "mandell_sims/input_thursday15_nions_beta1.out.nc",
    # #               "mandell_sims/input_thursday15_nions_beta1_light_electrons.out.nc",
    #                "mandell_sims/input_thursday16_nions_beta1_upar1.out.nc",
    # #               "mandell_sims/input_thursday13_smions_double_me.out.nc",
    # #               "mandell_sims/input_thursday14_smions_half_mi.out.nc",
    # #               #"mandell_sims/input_thursday11_smions_2.out.nc",
    # #               #"mandell_sims/input_thursday10_shions.out.nc",
    #                      ],
    #               [
    # #               "8",
    # #              "12",
    # #              "13",
    # #              "14",
    #               "15",
    # #              "15, light electrons",
    #               "16",
    #               ])

    # find_ksaw_properties_from_outnc(mandell_sf10_kperp001_new + ".out.nc")
    # examine_gs2_sim()

    # benchmark_stella_vs_gs2()
    benchmark_stella_vs_gs2_fapar1_fbpar0()
    # view_ncdf_variables_with_xarray("sims/stella_ky1_beta1_zero_gradients/input_fully_explicit.out.nc")
    # make_plots_for_movies("sims/stella_ky1_beta1_zero_gradients/input_fully_explicit.out.nc",
    #                        "stella", "movies/stella_ky1_beta1_zero_gradients")
    # make_plots_for_movies("sims/stella_ky1_beta1_zero_gradients_fapar0_fbpar0/input_fully_explicit_short.out.nc",
    #                        "stella", "movies/stella_ky1_beta1_zero_gradients_fapar0_fbpar0")
    # make_plots_for_movies("sims/gs2_ky1_beta1_zero_gradients_fapar1_fbpar1/input_long.out.nc",
    #                        "gs2", "movies/gs2_ky1_beta1_zero_gradients_fapar1_fbpar1")
