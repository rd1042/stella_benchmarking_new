"""Make plots for electromagnetic slab sims"""

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
# import ast

IMAGE_DIR = "./images/"

basic_em_sim = "stella_sims/input_slab_ky0.1_explicit"
mandell_beta1_kperp1 = "stella_sims/input_slab_ky0.1_explicit_mandell1"
mandell_beta1_kperp001_new = "mandell_sims/input_slab_ky0.01_beta1_new"
mandell_beta1_kperp1_new = "mandell_sims/input_slab_ky1_beta1_new"
mandell_sf1_kperp1_new = "mandell_sims/input_slab_ky1_sf1_new"
mandell_sf1_kperp001_new = "mandell_sims/input_slab_ky0.01_sf1_new"
mandell_sf0002_kperp001_new = "mandell_sims/input_slab_ky0.01_sf0.002_new"
mandell_sf00002_kperp001_new = "mandell_sims/input_slab_ky0.01_sf0.0002_new"
mandell_sf000002_kperp001_new = "mandell_sims/input_slab_ky0.01_sf0.00002_new"
mandell_sf10_kperp001_new = "mandell_sims/input_slab_ky0.01_sf10_new"
mandell_beta1_kperp1_long_t = "stella_sims/input_slab_ky0.1_explicit_mandell2"
mandell_beta1_kperp1_long_t_marconi = "mandell_sims/input_slab_ky0.1_beta1"

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

def compare_sims(outnc_longnames, labels):
    """ """

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    for idx, outnc_longname in enumerate(outnc_longnames):
        t, z, phi_vs_t, beta = extract_data_from_ncdf(outnc_longname, "t", 'zed', 'phi_vs_t', 'beta')
        phi_vs_t = phi_vs_t[:,0,:,0,0]
        print("len(t) = ", len(t))
        print("len(z) = ", len(z))
        print("phi_vs_t.shape = ", phi_vs_t.shape)
        # print("phi_vs_t = ", phi_vs_t)
        # print("phi_vs_t.imag = ", phi_vs_t.imag)
        # sys.exit()
        ax1.plot(t, phi_vs_t[:,int(len(z)*0.5)].real, label=(labels[idx] + "z=" + str(z[int(len(z)*0.5)]) ))
        ax2.plot(t, phi_vs_t[:,int(len(z)*0.5)].imag)
    ax2.set_xlabel("t")
    ax2.set_ylabel("Im(phi)")
    ax1.set_ylabel("Re(phi)")
    ax1.legend(loc="best")
    for ax in [ax1, ax2]:
        ax.grid(True)
    plt.show()

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

def benchmark_stella_vs_mandell():
    """ """
    phi_vs_t_longnames = glob.glob("mandell_sims/*.pickle")

    for file_longname in phi_vs_t_longnames:
        print("file_longname = ", file_longname)
        find_ksaw_properties_from_pickle(file_longname)
    return

def benchmark_stella_vs_mandell2():
    """ """
    phi_vs_t_longnames = glob.glob("mandell_sims_supermassive_ions/*.pickle")

    for file_longname in phi_vs_t_longnames:
        print("file_longname = ", file_longname)
        find_ksaw_properties_from_pickle(file_longname)
    return

def benchmark_stella_vs_gs2():
    """Compare similar simulations between stella and GS2"""

    stella_sim = "sims/stella_ky1_beta1_zero_gradients/input.out.nc"
    gs2_sim = "sims/gs2_ky1_beta1_zero_gradients/input.out.nc"
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
    compare_sims_with_apar([
                    "mandell_sims/input_thursday17_nions_beta1_write_apar.out.nc",
                    "mandell_sims/input_thursday18_wapar_upar.out.nc",
                    "mandell_sims/input_thursday19_2upar_beta1.out.nc",
                    "mandell_sims/input_thursday20_no_upar_beta10.out.nc",
                    "mandell_sims/input_thursday21_1upar_beta10.out.nc",
                    "mandell_sims/input_thursday21_2upar_beta10.out.nc",
                    "mandell_sims/input_thursday23_2upar_beta1_ky0.02.out.nc",
                    "mandell_sims/input_thursday23_2upar_beta1_ky0.1.out.nc",
                    "mandell_sims/input_thursday24_2upar_beta0.01_ky0.1.out.nc",
                    "mandell_sims/input_thursday22_2upar_beta1_ky1.out.nc",
                         ],
                  [
                  "17",
                  "18",
                  "19",
                  "20",
                  "21 1",
                  "21 2",
                  "23 1",
                  "23 2",
                  "24",
                  "22",
                  ])
    # find_ksaw_properties_from_outnc(mandell_sf10_kperp001_new + ".out.nc")
    examine_gs2_sim()
