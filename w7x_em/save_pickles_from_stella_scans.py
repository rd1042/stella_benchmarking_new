""" """

import numpy as np
import pickle
import glob
import sys

sys.path.append("../postprocessing_tools")
import helper_linear_sims as help_lin
from helper_ncdf_new import view_ncdf_variables, extract_data_from_ncdf_with_xarray
from extract_sim_data import get_omega_data, get_phiz_data, get_aparz_data, get_bparz_data
from plotting_helper import plot_phi_z_for_sim, plot_apar_z_for_sim, plot_bpar_z_for_sim


# sys.path.append("../postprocessing_tools_from_stella_dev")
# from extract_sim_data import get_omega_data_stella_from_ncdf
# from helper_ncdf_new import extract_data_from_ncdf_with_xarray
import re
import shutil
import os
import matplotlib.pyplot as plt

def get_phiz_data_stella(sim_longname):
    """ """
    final_fields_filename = sim_longname + ".final_fields"
    final_fields_file = open(final_fields_filename, "r")
    final_fields_data=np.loadtxt(final_fields_filename,dtype='float')
    final_fields_file.close()

    ## final_fields_data = z, z-zed0, aky, akx, real(phi), imag(phi), real(apar), imag(apar), z_eqarc-zed0, kperp2
    # Usually we're just looking at one mode; check how many unique kx and ky we have

    z = final_fields_data[:,0]

    aky = final_fields_data[:,2]; akx = final_fields_data[:,3]
    unique_aky = set(aky); unique_akx = set(akx)
    if len(unique_aky) > 1 or len(unique_akx) > 1:
        print("len(unique_aky), len(unique_akx) = ", len(unique_aky), len(unique_akx))
        print("Not currently supported")
        sys.exit()

    real_phi = final_fields_data[:,4]; imag_phi = final_fields_data[:,5]
    return z, real_phi, imag_phi

def get_fprim_tprim_ky(outnc_longname):
    """ """
    fprim, tprim, ky = extract_data_from_ncdf_with_xarray(outnc_longname, "fprim", "tprim", "ky")
    # print("fprim = ", float(fprim[0]))
    # print("tprim = ", float(tprim[0]))
    # print("ky = ", float(ky))
    # sys.exit()
    return float(fprim[0]), float(tprim[0]), float(ky)

def get_beta_from_outnc_longname(outnc_longname):
    """ """
    beta = extract_data_from_ncdf_with_xarray(outnc_longname, "beta")
    print("beta = ", beta)
    sys.exit()
    return float(beta)
    # return float(fprim[0]), float(tprim[0]), float(ky)

def construct_longlists_for_stella_fprim_tprim_ky_scan(folder_name, has_subfolders=True):
    """For a STEP folder, extract psin, ky or n, and some other variable
    (currently intended to be either kx or theta0). Put this data
    into long flat lists, and create corresponding lists for
    frequency and growth rate. Also create plots:
    (1) Converence plots
    (2) field plots

    If with_epar=True, also calculate epar, and include epar data in the field plots."""

    # Find all the simulation files - look for everything with a .in
    # infile_longnames = glob.glob(folder_name + "/*/*.in")
    # Alas, gs2's master now has a "<name>.used_inputs.in", so the above is a
    # bad idea. Instead, look for .out.nc files.
    if has_subfolders:
        subfolder_longnames = glob.glob(folder_name + "/*/")
    else:
        outnc_longnames = glob.glob(folder_name + "/*.out.nc")

    fprim_longlist = []
    tprim_longlist = []
    ky_longlist = []

    omega_file_growth_rate_longlist = []
    safe_growth_rate_longlist = []
    freq_longlist = []
    bad_sims_list = []  # Records sims which somehow failed
    unconverged_sim_list = []
    if has_subfolders:
        for subfolder_longame in subfolder_longnames:
            # Find the .out.nc files
            outnc_longnames = glob.glob(subfolder_longame+ "*.out.nc")
            if len(outnc_longnames) > 1:
                # Get the latest sim - this will be input_X , where X=len(outnc_longnames)-1
                outnc_longname = subfolder_longame + "input_" + str(len(outnc_longnames)-1) + ".out.nc"
            else:
                outnc_longname = outnc_longnames[0]

            # If has_subfolders, Find the "latest" outnc_longname in the folder
            # Get the param vals and append to their lists
            fprim_val, tprim_val, ky_val = get_fprim_tprim_ky(outnc_longname)

            fprim_longlist.append(fprim_val)
            tprim_longlist.append(tprim_val)
            ky_longlist.append(ky_val)

            freqom_final, gammaom_final, gamma_safe, sim_converged = get_omega_from_subfolder_and_make_plot(outnc_longname)
            if not sim_converged:
                unconverged_sim_list.append(outnc_longname)
            omega_file_growth_rate_longlist.append(gammaom_final)
            safe_growth_rate_longlist.append(gamma_safe)
            freq_longlist.append(freqom_final)
    else:
        print("Not has_subfolders not supported, quitting")
        sys.exit()

    print("len(unconverged_sim_list) = ", len(unconverged_sim_list))
    print("unconverged_sim_list = ", unconverged_sim_list)

    return (fprim_longlist, tprim_longlist, ky_longlist,
            omega_file_growth_rate_longlist, safe_growth_rate_longlist,
            freq_longlist)

def construct_longlists_for_stella_beta_scan(folder_name, has_subfolders=True):
    """Also create plots:
    (1) Converence plots
    (2) field plots

    If with_epar=True, also calculate epar, and include epar data in the field plots."""

    # Find all the simulation files - look for everything with a .in
    # infile_longnames = glob.glob(folder_name + "/*/*.in")
    # Alas, gs2's master now has a "<name>.used_inputs.in", so the above is a
    # bad idea. Instead, look for .out.nc files.
    if has_subfolders:
        subfolder_longnames = glob.glob(folder_name + "/*/")
    else:
        outnc_longnames = glob.glob(folder_name + "/*.out.nc")

    beta_longlist = []

    omega_file_growth_rate_longlist = []
    safe_growth_rate_longlist = []
    freq_longlist = []
    bad_sims_list = []  # Records sims which somehow failed
    unconverged_sim_list = []
    if has_subfolders:
        for subfolder_longame in subfolder_longnames:
            # Find the .out.nc files
            outnc_longnames = glob.glob(subfolder_longame+ "*.out.nc")
            if len(outnc_longnames) > 1:
                # Get the latest sim - this will be input_X , where X=len(outnc_longnames)-1
                outnc_longname = subfolder_longame + "input_" + str(len(outnc_longnames)-1) + ".out.nc"
            else:
                outnc_longname = outnc_longnames[0]

            # If has_subfolders, Find the "latest" outnc_longname in the folder
            # Get the param vals and append to their lists
            beta_val = get_beta_from_outnc_longname(outnc_longname)

            beta_longlist.append(beta_val)

            freqom_final, gammaom_final, gamma_safe, sim_converged = get_omega_from_subfolder_and_make_plot(outnc_longname)
            if not sim_converged:
                unconverged_sim_list.append(outnc_longname)
            omega_file_growth_rate_longlist.append(gammaom_final)
            safe_growth_rate_longlist.append(gamma_safe)
            freq_longlist.append(freqom_final)
    else:
        print("Not has_subfolders not supported, quitting")
        sys.exit()

    print("len(unconverged_sim_list) = ", len(unconverged_sim_list))
    print("unconverged_sim_list = ", unconverged_sim_list)

    return (beta_longlist,
            omega_file_growth_rate_longlist, safe_growth_rate_longlist,
            freq_longlist)

def make_omega_arrays_from_longlists(fprim_longlist, tprim_longlist, ky_longlist,
                    omega_file_growth_rate_longlist, safe_growth_rate_longlist,
                    freq_longlist):
    """ """
    # Construct the array
    unique_fprim = sorted(set(fprim_longlist))
    unique_tprim = sorted(set(tprim_longlist))
    unique_ky = sorted(set(ky_longlist))
    omega_file_gamma_fprim_tprim_ky_array = np.zeros((len(unique_fprim), len(unique_tprim),
                                    len(unique_ky)))
    safe_gamma_fprim_tprim_ky_array = np.zeros((len(unique_fprim), len(unique_tprim),
                                    len(unique_ky)))
    omega_fprim_tprim_ky_array = np.zeros((len(unique_fprim), len(unique_tprim),
                                    len(unique_ky)))
    # Initialise with NaNs:
    for my_array in [omega_file_gamma_fprim_tprim_ky_array,
                     safe_gamma_fprim_tprim_ky_array,
                     omega_fprim_tprim_ky_array]:
        my_array[:,:,:] = np.NaN

    # Populate the array
    for i in range(0, len(fprim_longlist)):
        fprim = fprim_longlist[i]
        tprim = tprim_longlist[i]
        ky = ky_longlist[i]
        omega_file_growth_rate = omega_file_growth_rate_longlist[i]
        safe_growth_rate = safe_growth_rate_longlist[i]
        freq = freq_longlist[i]

        fprim_idx = unique_fprim.index(fprim)
        tprim_idx = unique_tprim.index(tprim)
        ky_idx = unique_ky.index(ky)
        omega_file_gamma_fprim_tprim_ky_array[fprim_idx, tprim_idx, ky_idx] = omega_file_growth_rate
        safe_gamma_fprim_tprim_ky_array[fprim_idx, tprim_idx, ky_idx] = safe_growth_rate
        omega_fprim_tprim_ky_array[fprim_idx, tprim_idx, ky_idx] = freq

    print("empty spaces = ", (len(unique_fprim) * len(unique_tprim) * len(unique_ky) - len(fprim_longlist)))

    return unique_fprim, unique_tprim, unique_ky, omega_file_gamma_fprim_tprim_ky_array, safe_gamma_fprim_tprim_ky_array, omega_fprim_tprim_ky_array

def postprocess_fprim_tprim_ky_scan(folder_name, pickle_string=None):
    """ """

    (fprim_longlist, tprim_longlist, ky_longlist,
     omega_file_growth_rate_longlist, safe_growth_rate_longlist,
            freq_longlist) = construct_longlists_for_stella_fprim_tprim_ky_scan(folder_name,
                                                has_subfolders=True)
    (unique_fprim, unique_tprim, unique_ky, omega_file_gamma_fprim_tprim_ky_array,
     safe_gamma_fprim_tprim_ky_array,
     omega_fprim_tprim_ky_array) = make_omega_arrays_from_longlists(fprim_longlist, tprim_longlist, ky_longlist,
                                            omega_file_growth_rate_longlist, safe_growth_rate_longlist,
                                            freq_longlist)
    print("unique_fprim = ", unique_fprim)
    print("unique_tprim = ", unique_tprim)
    print("unique_ky = ", unique_ky)
    print("omega_file_gamma_fprim_tprim_ky_array = ", omega_file_gamma_fprim_tprim_ky_array)
    print("safe_gamma_fprim_tprim_ky_array = ", safe_gamma_fprim_tprim_ky_array)
    print("omega_fprim_tprim_ky_array = ", omega_fprim_tprim_ky_array)
    # Write to pickle
    pickle_shortname = "omega_fprim_tprim_ky_array.pickle"
    pickle_longname = folder_name + "/" + pickle_shortname
    if pickle_string is None:
        pickle_string = folder_name
    pickle_file = open(pickle_longname, "wb")
    pickle.dump([pickle_string, unique_fprim, unique_tprim, unique_ky,
        omega_file_gamma_fprim_tprim_ky_array, safe_gamma_fprim_tprim_ky_array,
        omega_fprim_tprim_ky_array], pickle_file)
    pickle_file.close()

    return

def postprocess_beta_scan(folder_name, pickle_string=None):
    """ """

    (beta_longlist,
     omega_file_growth_rate_longlist, safe_growth_rate_longlist,
            freq_longlist) = construct_longlists_for_stella_beta_scan(folder_name,
                                                has_subfolders=True)
    beta_array = np.array(beta_longlist)
    omega_file_gamma_fprim_tprim_ky_array = np.array(omega_file_gamma_beta_array)
    safe_gamma_fprim_tprim_ky_array = np.array(safe_gamma_beta_array)
    omega_fprim_tprim_ky_array = np.array(omega_beta_array)
    print("beta_array = ", beta_array)
    print("omega_file_gamma_fprim_tprim_ky_array = ", omega_file_gamma_beta_array)
    print("safe_gamma_fprim_tprim_ky_array = ", safe_gamma_beta_array)
    print("omega_fprim_tprim_ky_array = ", omega_beta_array)
    # Write to pickle
    pickle_shortname = "omega_beta_array.pickle"
    pickle_longname = folder_name + "/" + pickle_shortname
    if pickle_string is None:
        pickle_string = folder_name
    pickle_file = open(pickle_longname, "wb")
    pickle.dump([pickle_string, beta_array,
        omega_file_gamma_fprim_tprim_ky_array, safe_gamma_fprim_tprim_ky_array,
        omega_fprim_tprim_ky_array], pickle_file)
    pickle_file.close()

    return

def get_omega_from_subfolder_and_make_plot(outnc_longname):
    """For a folder, get Omega and beta, and save to pickle"""

    sim_shortname = re.split("/", outnc_longname)[-1]
    sim_shortname = re.split(".out.nc", sim_shortname)[0]
    sim_longname = re.split(".out.nc", outnc_longname)[0]
    # view_ncdf_variables(outnc_longname)
    time, freqom_final, gammaom_final, freqom, gammaom, gamma_stable = get_omega_data(sim_longname, "stella")
    # These might not be finite - if not, from the finite versions
    # Create a flag to see if finite
    if np.isfinite(freqom).all():
        freqom_all_finite = True
    else:
        freqom_all_finite = False
    if np.isfinite(gammaom).all():
        gammaom_all_finite = True
    else:
        gammaom_all_finite = False
    if np.isfinite(gamma_stable).all():
        gamma_stable_all_finite = True
    else:
        gamma_stable_all_finite = False
    freqom_finite_idxs = np.isfinite(freqom)
    freqom_finite = freqom[freqom_finite_idxs]
    time_freqom = time[freqom_finite_idxs]
    gammaom_finite_idxs = np.isfinite(gammaom)
    gammaom_finite = gammaom[gammaom_finite_idxs]
    time_gammaom = time[gammaom_finite_idxs]
    gamma_stable_finite_idxs = np.isfinite(gamma_stable)
    gamma_stable_finite = gamma_stable[gamma_stable_finite_idxs]

    freqom_tolerance = 1E-2
    freqom_frac_tolerance = 2E-2
    gammaom_tolerance = 1E-2
    gammaom_frac_tolerance = 2E-2
    gamma_stable_tolerance = 1E-2
    gamma_frac_stable_tolerance = 2E-2
    time_tolerance = 1E-5
    sim_converged = True
    ## Check if freqom, gammaom, gamma_stable are converged
    if ((np.abs(freqom_finite[-1] - freqom_finite[-2]) > freqom_tolerance) and
            (np.abs((freqom_finite[-1] - freqom_finite[-2])/freqom_finite[-1]) > freqom_frac_tolerance)):
        sim_converged = False
    if ((np.abs(gammaom_finite[-1] - gammaom_finite[-2]) > gammaom_tolerance) and
            (np.abs((gammaom_finite[-1] - gammaom_finite[-2])/gammaom_finite[-1]) > gammaom_frac_tolerance)):
        sim_converged = False
    if ((np.abs(gamma_stable_finite[-1] - gamma_stable_finite[-2]) > gamma_stable_tolerance) and
            (np.abs((gamma_stable_finite[-1] - gamma_stable_finite[-2])/gamma_stable_finite[-1]) > gamma_frac_stable_tolerance)):
        sim_converged = False
    if ((np.abs(gamma_stable_finite[-1] - gammaom_finite[-1]) > gamma_stable_tolerance) and
            (np.abs((gamma_stable_finite[-1] - gammaom_finite[-1])/gamma_stable_finite[-1]) > gamma_frac_stable_tolerance)):
        sim_converged = False
    if not sim_converged:
        # Either run for longer, or drop the timestep
        if (freqom_all_finite) and (gammaom_all_finite) and (gamma_stable_all_finite):
            # Just need to increase time
            increase_sim_time(sim_longname)
        else:
            # If time[-1] is the same for finite_freqom etc., then apparently the non-finite values
            # turn up somewhere else
            if ( (abs((time[-1] - time_freqom[-1])) < time_tolerance) and
                 (abs((time[-1] - time_gammaom[-1])) < time_tolerance) and
                 (abs((time[-1] - time_gamma_stable[-1])) < time_tolerance) ):
                print("Don't know why the sim has failed to converged")
                print("sim_shortname = ", sim_shortname)
                print("(abs((time[-1] - time_freqom[-1])), " +
                      "(abs((time[-1] - time_gammaom[-1])), " +
                      "(abs((time[-1] - time_gamma_stable[-1])) = ",
                      abs(time[-1] - time_freqom[-1]), abs(time[-1] - time_gammaom[-1]),
                      abs(time[-1] - time_gamma_stable[-1]) )
                print("freqom = ", freqom)
                print("gammaom = ", gammaom)
                print("gamma_stable = ", gamma_stable)
                sys.exit()
            else:
                decrease_dt(sim_longname)
    #try:
    z, real_phi, imag_phi = get_phiz_data_stella(sim_longname)
    if (np.isfinite(real_phi)).all() and (np.isfinite(imag_phi)).all():
        abs_phi = np.abs(real_phi + 1j*imag_phi)
    else:
        abs_phi = np.ones(len(z)) * np.nan


    # make_omega_time_plot_for_stella(sim_longname, time, freqom, gammaom, gamma_stable)
    make_phi_z_plot_for_stella(sim_longname, z, abs_phi)

    return freqom[-1], gammaom[-1], gamma_stable[-1], sim_converged

def make_phi_z_plot_for_stella(sim_longname, z, abs_phi):
    """ """
    fig = plt.figure(figsize=(12,12))
    ax1 = fig.add_subplot(111)
    if (np.isfinite(abs_phi)).all():
        ax1.plot(z/np.pi, abs_phi/np.max(abs_phi), c="black")
    ax1.grid(True)
    ax1.set_xlabel(r"$z/\pi$")
    ax1.set_ylabel(r"$\vert \phi \vert$")
    plt.tight_layout()
    # Work out the save name
    # sim_longname looks like sim_longname=<parent_folder>/<parameters>/input(_2,_3,...)
    parameter_simname = re.split("/", sim_longname)[-2]
    parent_folder = re.split("/"+parameter_simname, sim_longname)[0]
    sim_shortname = re.split("/", sim_longname)[-1]

    save_name = parent_folder + "/" + parameter_simname
    if "_" in sim_shortname:
        numstr = re.split("_", sim_shortname)[-1]
        save_name = save_name + "_"
    save_name = save_name + ".png"
    plt.savefig(save_name)
    plt.close()


    return

def increase_sim_time(sim_longname):
    """Make a copy of the sim folder with an appropriate suffix, with the same
       input file and runscript but different nstep"""
    sim_folder_name = re.split("/input", sim_longname)[0]
    sim_shortname = re.split("/", sim_folder_name)[-1]
    folder_prefix = re.split(sim_shortname, sim_folder_name)[0]
    simfile_shortname = re.split(sim_folder_name+"/", sim_longname)[-1]
    ## Find out if the simfile is "input", "input_1" etc.
    split_simfile = re.split("_", simfile_shortname)
    if len(split_simfile) == 1:
        # it's "input"
        input_counter = 1
    elif len(split_simfile) ==2:
        input_counter = int(split_simfile[-1])
        input_counter += 1

    new_simfile_shortname = "input_" + str(input_counter)
    new_runscript_longname = sim_folder_name + "/run_stella_" + str(input_counter) + ".sh"
    template_runscript = open(sim_folder_name + "/run_stella.sh", "r")
    runscript_text = template_runscript.read()
    template_runscript.close()
    search_string = "input.in"
    replace_string = new_simfile_shortname + ".in"
    runscript_text = re.sub(search_string, replace_string, runscript_text)
    new_runscript = open(new_runscript_longname, "w")
    new_runscript.write(runscript_text)
    new_runscript.close()
    #shutil.copy(sim_folder_name + "/run_stella.sh", new_runscript_longname)
    ### Copy over the input file, but increasing nstep by 50%

    template_file = open(sim_longname+".in", 'r+')
    filetext = template_file.read()
    template_file.close()
    lines = re.split("\n", filetext)
    for line_idx, line in enumerate(lines):

        if "nstep" in line:
            nstep_old = float((re.split("nstep = ", line)[-1]).strip())
    new_nstep_str = str(int(nstep_old*1.8))
    search_string = "nstep =.*"
    replace_string = "nstep = " + new_nstep_str
    #print("search_string = ", search_string)
    #print("replace_string = ", replace_string)
    filetext = re.sub(search_string, replace_string, filetext)

    new_input_file = sim_folder_name + "/input_" + str(input_counter) + ".in"
    new_infile = open(new_input_file, "w")
    new_infile.write(filetext)
    new_infile.close()

    return

def decrease_dt(sim_longname):
    """Make a copy of the sim folder with an appropriate suffix, with the same
       input file and runscript but different nstep"""
    sim_folder_name = re.split("/input", sim_longname)[0]
    sim_shortname = re.split("/", sim_folder_name)[-1]
    folder_prefix = re.split(sim_shortname, sim_folder_name)[0]
    simfile_shortname = re.split(sim_folder_name, sim_longname)[-1]
    ## Find out if the simfile is "input", "input_1" etc.
    split_simfile = re.split("_", simfile_shortname)
    if len(split_simfile) == 1:
        # it's "input"
        input_counter = 1
    elif len(split_simfile) ==2:
        input_counter = float(split_simfile[-1])
        input_counter += 1
    print("sim_longname = ", sim_longname)
    new_simfile_shortname = "input_" + str(input_counter)
    new_runscript_longname = sim_folder_name + "/run_stella_" + str(input_counter) + ".sh"
    template_runscript = open(sim_folder_name + "/run_stella.sh", "r")
    runscript_text = template_runscript.read()
    template_runscript.close()
    search_string = simfile_shortname + ".in"
    replace_string = new_simfile_shortname + ".in"
    runscript_text = re.sub(search_string, replace_string, runscript_text)
    new_runscript = open(new_runscript_longname, "w")
    new_runscript.write(runscript_text)
    new_runscript.close()

    ### Copy over the input file, but reducing code_dt by a factor of 2

    template_file = open(sim_longname+".in", 'r+')
    filetext = template_file.read()
    template_file.close()
    lines = re.split("\n", filetext)
    for line_idx, line in enumerate(lines):

        if "delt =" in line:
            code_dt = float((re.split("delt = ", line)[-1]).strip())
    new_code_dt_str = str(code_dt/2)
    search_string = "delt =.*"
    replace_string = "delt = " + new_code_dt_str
    #print("search_string = ", search_string)
    #print("replace_string = ", replace_string)
    filetext = re.sub(search_string, replace_string, filetext)

    new_input_file = sim_folder_name + "/input_" + str(input_counter) + ".in"
    new_infile = open(new_input_file, "w")
    new_infile.write(filetext)
    new_infile.close()
    return
