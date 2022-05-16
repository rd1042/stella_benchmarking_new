"""Script to read netCDF output from the linear cyclone ITG base case and fit the
growth rate

Author: Bob Davies"""

import helper_ncdf_new as help_ncdf
import glob
import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
import re
import math
import shutil
import sys
import matplotlib as mpl
import ncdf2dict

def get_param_from_input_file(infile_longname, param_name):
    """Open an input file, search for param_name, and return the value"""
    myfile = open(infile_longname, "r")
    filetext = myfile.read()
    myfile.close()
    param_lines = []
    param_val_strs = []
    for line in re.split("\n", filetext):
        if param_name in line:
            param_lines.append(line)
            param_val_strs.append((re.split("=", line)[-1]).strip())

    if len(param_val_strs) > 1:
        print("param_lines = ", param_lines)
        print("param_val_strs = ", param_val_strs)
        return param_val_strs
    else:
        return param_val_strs[0]

def plot_param_scan(folder_name):
    """Plot the parameter scan"""
    print("folder_name = ", folder_name)
    (param_name, param_vals, frequency_list, freq_errors, growth_rate_list,
     growth_rate_errors) = extract_params_from_pickle(folder_name)
    # Make the plot
    make_param_scan_plot(param_name, param_vals, frequency_list, freq_errors,
                         growth_rate_list, growth_rate_errors)

    return

def plot_growth_rate_for_param_scan(folder_name):
    """Plot the parameter scan"""
    print("folder_name = ", folder_name)
    (param_name, param_vals, growth_rate_list, growth_rate_errors) = extract_growth_rate_from_pickle(folder_name)
    # Make the plot
    make_growth_rate_plot(param_name, param_vals, growth_rate_list, growth_rate_errors, folder_name)

    return

def make_growth_rate_plot(param_name, param_list, growth_rate_list, growth_rate_errors, sim_name, nunits="normalised GS2 units"):

    # Plot the figure
    fig = plt.figure()
    ax2 = fig.add_subplot(111)
    ax2.scatter(param_list, growth_rate_list, c='black', marker='x')
    ax2.errorbar(param_list, growth_rate_list, yerr=growth_rate_errors, fmt='none',
    c='r', capsize=3)
    print("[min(param_list)[0], max(param_list)[0]] = ", [min(param_list)[0][0], max(param_list)[0][0]])
    ax2.plot([min(param_list)[0][0], max(param_list)[0][0]], [0, 0], c="black")
    ax2.set_xlabel(param_name + " parameter")
    ax2.set_ylabel("Growth rate (" + nunits + ")")
    ax2.set_xlim(min(param_list)[0][0], max(param_list)[0][0])
    plt.title(sim_name)
    plt.show()
    return

def make_param_scan_plot(param_name, param_list, frequency_list, freq_errors,
                         growth_rate_list, growth_rate_errors, nunits="normalised GS2 units"):

    # Plot the figure
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.scatter(param_list, frequency_list, c='black', marker='x')
    ax1.errorbar(param_list, frequency_list, yerr=freq_errors, fmt='none',
    c='r', capsize=3)
    ax1.set_xlabel(param_name +  " parameter")
    ax1.set_ylabel("Frequency (" + nunits + ")")

    ax2 = fig.add_subplot(212)
    ax2.scatter(param_list, growth_rate_list, c='black', marker='x')
    ax2.errorbar(param_list, growth_rate_list, yerr=growth_rate_errors, fmt='none',
    c='r', capsize=3)
    ax2.set_xlabel(param_name + " parameter")
    ax2.set_ylabel("Growth rate (" + nunits + ")")
    plt.show()
    return

def calculate_omega_for_param_scan(folder_name, param, param_key, file_regex=False,
                                   path=None, use_lin_filepath=False, **kwargs):
    """Extract output files from selected folder. For each output file, find the
    value of the scanned parameter and calculate Omega."""
    print("folder_name = ", folder_name)
    # print("Check your slashes")
    # sys.exit()
    # Extract the files
    if path is None:
        if use_lin_filepath:
            path = LIN_SIM_PATH + folder_name + "/"
        elif file_regex == False:
            path = folder_name + '/'
        else:
            path = folder_name + '/' + file_regex
    output_file_list = glob.glob(path + "*.out.nc")
    param_list = []
    frequency_list = []
    growth_rate_list = []
    freq_errors = []
    growth_rate_errors = []
    if len(output_file_list) == 0:
        print("output_file_list = ", output_file_list)
        sys.exit("Empty output file list, aborting")
    for i, output_file in enumerate(output_file_list):

        # Extract the name of the simulation from its full path
        m = re.search(path +'(.+?)\.out\.nc', output_file)
        found = m.group(1)
        save_loc = re.split((found+".out.nc"), output_file)[0]
        print("save_loc = ", save_loc)
        fitted_growth_rate, growth_rate_error, converged_frequency, freq_error = calculate_omega_and_plot_for_single(output_file,
                        figtitle=found, folder_name=folder_name, save_path=save_loc, **kwargs)
        frequency_list.append(converged_frequency)
        freq_errors.append(freq_error)
        growth_rate_list.append(fitted_growth_rate)
        growth_rate_errors.append(growth_rate_error)
        print("output_file, fitted_growth_rate = ", output_file, fitted_growth_rate)
        try:    # Bob: should fix this
            param_list.append(help_ncdf.extract_data_from_ncdf(output_file, param_key))
        except KeyError:
            print("Error!")
            sys.exit()
    return (param_list, frequency_list, growth_rate_list,
            freq_errors, growth_rate_errors)

def calculate_omega_and_plot_for_single(outnc_longname, plot_growth_rates=False, save_path=None,
                        folder_name=".", figtitle="Fig title"):
    """ """

    if (save_path is None) and (plot_growth_rates):
        print("Warning! No save path for plots specified")
        save_path = "./"

    # Plot and save the growth rates
    if plot_growth_rates:
        help_ncdf.plot_growth_rates(outnc_longname, figtitle, save_loc=save_path)
    # Find the value of the parameter being scanned

    # Find the frequency and growth rate, and append to the lists
    return help_ncdf.calculate_omega(outnc_longname)

def calculate_omega_for_2d_param_scan(folder_name, param1, param2, param_key1, param_key2,
                                        use_lin_filepath=True, plot_growth_rates=True):
    """Extract output files from selected folder. For each output file, find the
    value of the scanned parameter and calculate Omega."""

    # Extract the files
    if use_lin_filepath:
        path = LIN_SIM_PATH + folder_name + '/'
    else:
        path = folder_name + '/'
    output_file_list = glob.glob(path + '*.out.nc')

    param1_list = []
    param2_list = []
    frequency_list = []
    growth_rate_list = []
    freq_errors = []
    growth_rate_errors = []
    if len(output_file_list) == 0:
        print("output_file_list = ", output_file_list)
        sys.exit("Empty output file list, aborting")
    for output_file in output_file_list:

        # Extract the name of the simulation from its full path
        m = re.search(path +'(.+?)\.out\.nc', output_file)
        found = m.group(1)
        print("found = ", found)

        # Plot and save the growth rates
        if plot_growth_rates:
            help_ncdf.plot_growth_rates(output_file, found, save_loc=path)
        # Find the value of the parameter being scanned
        # print("param_key1 = ", param_key1)
        # print("param1 = ", help_ncdf.extract_data_from_ncdf(output_file, param_key1))
        # sys.exit()
        param1_list.append(help_ncdf.extract_data_from_ncdf(output_file, param_key1))
        param2_list.append(help_ncdf.extract_data_from_ncdf(output_file, param_key2))
        # Find the frequency and growth rate, and append to the lists
        fitted_growth_rate, growth_rate_error, converged_frequency, freq_error= (
        help_ncdf.calculate_omega(output_file))
        frequency_list.append(converged_frequency)
        freq_errors.append(freq_error)
        growth_rate_list.append(fitted_growth_rate)
        growth_rate_errors.append(growth_rate_error)
        print("output_file, fitted_growth_rate = ", output_file, fitted_growth_rate)

    return (param1_list, param2_list, frequency_list, growth_rate_list,
            freq_errors, growth_rate_errors)

def calculate_growth_rate_for_param_scan(folder_name, param, param_key):
    """Extract output files from selected folder. For each output file, find the
    value of the scanned parameter and calculate Omega."""

    # Extract the files
    path = LIN_SIM_PATH + folder_name + '/'
    output_file_list = glob.glob(path + '*.out.nc')

    param_list = []
    growth_rate_list = []
    growth_rate_errors = []

    for output_file in output_file_list:

        # Extract the name of the simulation from its full path
        m = re.search('\.\./\.\./gs2_sims_linear/' + folder_name +'/(.+?)\.out\.nc', output_file)
        found = m.group(1)

        # Plot and save the growth rates
        #help_ncdf.plot_growth_rates(output_file, found, save_loc=path)
        # Find the value of the parameter being scanned
        param_list.append(help_ncdf.extract_data_from_ncdf(output_file, param_key))
        # Find the frequency and growth rate, and append to the lists
        fitted_growth_rate, growth_rate_error = (help_ncdf.lazy_calculate_omega(output_file))
        growth_rate_list.append(fitted_growth_rate)
        growth_rate_errors.append(growth_rate_error)

    return (param_list, growth_rate_list, growth_rate_errors)

def save_growth_rate_for_param_scan(folder_name, param, param_key):
    """Calculate growth rate and frequency for each simulation in a given parameter
    scan folder, and save the data to a file."""

    (param_list, growth_rate_list, growth_rate_errors) = calculate_growth_rate_for_param_scan(folder_name, param, param_key)

    # Create dictionary of values
    omega_dict = {"parameter name" : param,
            param : param_list,
            "growth rate" : growth_rate_list,
            "growth rate error" : growth_rate_errors}

    filename = LIN_SIM_PATH + folder_name + '/growth_rate_values.pickle'

    # Check if file already exists. If so, delete
    if os.path.isfile(filename):
        print("omega values file already exists, deleting")
        os.system('rm ' + filename)

    myfile = open(filename, 'wb+')
    pickle.dump(omega_dict, myfile)
    myfile.close()
    return

def save_omega_for_param_scan(folder_name, param1, param_key1, param2=False, param_key2=False,
                              file_regex=False, use_lin_filepath=True, plot_growth_rates=True):
    """Calculate growth rate and frequency for each simulation in a given parameter
    scan folder, and save the data to a file."""

    if param2 == False:
        (param_list, frequency_list, growth_rate_list,
                freq_errors, growth_rate_errors) = calculate_omega_for_param_scan(folder_name, param1, param_key1,
                                                    file_regex, use_lin_filepath=use_lin_filepath,
                                                    plot_growth_rates=plot_growth_rates)

        # Create dictionary of values
        omega_dict = {"parameter name" : param1,
                param1 : param_list,
                "growth rate" : growth_rate_list,
                "frequency" : frequency_list,
                "growth rate error" : growth_rate_errors,
                "frequency error" : freq_errors}

    else:
        (param1_list, param2_list, frequency_list, growth_rate_list,
        freq_errors, growth_rate_errors) = calculate_omega_for_2d_param_scan(folder_name, param1, param2, param_key1, param_key2,
                                                use_lin_filepath=use_lin_filepath, plot_growth_rates=plot_growth_rates)

        # Create dictionary of values
        omega_dict = {"parameter names" : [param1, param2],
        param1 : param1_list,
        param2 : param2_list,
        "growth rate" : growth_rate_list,
        "frequency" : frequency_list,
        "growth rate error" : growth_rate_errors,
        "frequency error" : freq_errors}

    if use_lin_filepath:
        filename = LIN_SIM_PATH + folder_name + '/omega_values.pickle'
    else:
        filename = folder_name + '/omega_values.pickle'

    # Check if file already exists. If so, delete
    if os.path.isfile(filename):
        print("omega values file already exists, deleting")
        os.system('rm ' + filename)

    myfile = open(filename, 'wb+')
    pickle.dump(omega_dict, myfile)
    myfile.close()
    return

def save_omega_alpha_for_salpha(folder_name):
    """ """

    # Find the folders in the parent directory
    child_folder_list = glob.glob(LIN_SIM_PATH + "/" + folder_name + "/*/")
    # For each folder, open an output file and find fprim, tprim, then save
    # omega(bprim, beta) for each.

    for child_folder_longname in child_folder_list:
        print("longname = " + child_folder_longname)
        # We need this to look like "folder_name/child_folder" - split based
        # on "/" to get the child folder, then construct the folder name
        child_folder_shortname = re.split("/", child_folder_longname)[-2]
        child_folder = folder_name + "/" + child_folder_shortname
        print("child_folder = ", child_folder)
        save_omega_for_param_scan(child_folder, 'beta', 'beta')
        #plot_omega_for_bprim_scans([child_folder], [child_folder], legend_font=16.0)
        #sys.exit()

    return

def load_omega_for_param_scan(folder_name):
    """Loads data from the pickled file of omega for a parameter scan"""

    print("folder_name = ", folder_name)
    filename = LIN_SIM_PATH + folder_name + '/omega_values.pickle'
    myfile = open(filename, 'rb')
    print("filename = ", filename)
    dict = pickle.load(myfile)

    return dict

def load_growth_rate_for_param_scan(folder_name):
    """Loads data from the pickled file of omega for a parameter scan"""

    print("folder_name = ", folder_name)
    filename = LIN_SIM_PATH + folder_name + '/growth_rate_values.pickle'
    myfile = open(filename, 'rb')
    dict = pickle.load(myfile)

    return dict

def plot_omega_for_beta_scans(folder_list, folder_labels, normalisations=[], legend_font=8, title=False, use_cols=False):

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    cols=["blue", "orange", "green", "black"]
    ax1.tick_params('x', which='both', labelsize=8, direction="in")
    ax1.grid(True, which="both", axis="both", linestyle="--")
    ax2.set_xlabel(r"beta")
    ax2.set_ylabel(r"$\gamma/(v_{th,ref}/L_{ref})$", rotation=0, labelpad=25)
    ax1.set_ylabel(r"$\omega/(v_{th,ref}/L_{ref})$", rotation=0, labelpad=25)
    ax2.grid(True, which="both", axis="both", linestyle="--")

    if title != False:
        fig.suptitle(title, fontsize=16)

    for i, folder in enumerate(folder_list):
        scan_dict = load_omega_for_param_scan(folder)

        scan_beta = scan_dict['beta']; scan_frequency=scan_dict['frequency']; scan_growth_rate = scan_dict['growth rate']
        print("scan_beta = ", scan_beta)
        scan_beta, scan_frequency, scan_growth_rate = list(zip(*sorted(zip(scan_beta, scan_frequency, scan_growth_rate))))
        if len(normalisations) > 0:
            scan_frequency = np.array(scan_frequency) * normalisations[i]
            scan_growth_rate = np.array(scan_growth_rate) * normalisations[i]

        if use_cols == True:
            ax1.plot(scan_beta, scan_frequency, label=folder_labels[i], c=cols[i], lw = 3.0)
        else:
            ax1.plot(scan_beta, scan_frequency, label=folder_labels[i], lw = 3.0)
        #ax1.set_xlabel("beta parameter")
        #ax1.set_ylabel("Frequency (normalised GS2 units)")
        ax1.legend(loc='best')

        if use_cols == True:
            ax2.plot(scan_beta, scan_growth_rate, label=folder_labels[i], c=cols[i], lw = 3.0)
        else:
            ax2.plot(scan_beta, scan_growth_rate, label=folder_labels[i], lw = 3.0)
        #ax2.set_xlabel("beta parameter")
        #ax2.set_ylabel("Growth rate (normalised GS2 nuits)")
        ax2.legend(loc='best', prop={'size': legend_font})
    plt.tight_layout()
    plt.show()

    return

def plot_omega_for_bprim_scans(folder_list, folder_labels, normalisations=[], legend_font=8, title=False, use_cols=False):
    """Massive cheat, we're just working out beta' because dens, temp, fprim, tprim are constant. Should
    fix this."""
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    cols=["blue", "orange", "green", "black"]
    ax1.tick_params('x', which='both', labelsize=8, direction="in")
    ax1.grid(True, which="both", axis="both", linestyle="--")
    ax2.set_xlabel(r"-d$\beta/$d$\rho$")
    ax2.set_ylabel(r"$\gamma/(v_{th,ref}/L_{ref})$", rotation=0, labelpad=25)
    ax1.set_ylabel(r"$\omega/(v_{th,ref}/L_{ref})$", rotation=0, labelpad=25)
    ax2.grid(True, which="both", axis="both", linestyle="--")

    if title != False:
        fig.suptitle(title, fontsize=16)

    for i, folder in enumerate(folder_list):
        scan_dict = load_omega_for_param_scan(folder)

        scan_beta = scan_dict['beta']; scan_frequency=scan_dict['frequency']; scan_growth_rate = scan_dict['growth rate']
        print("scan_beta = ", scan_beta)
        scan_beta, scan_frequency, scan_growth_rate = list(zip(*sorted(zip(scan_beta, scan_frequency, scan_growth_rate))))

        scan_bprim = np.array(scan_beta) * 18.2 ## FIX!
        if len(normalisations) > 0:
            scan_frequency = np.array(scan_frequency) * normalisations[i]
            scan_growth_rate = np.array(scan_growth_rate) * normalisations[i]

        if use_cols == True:
            ax1.plot(scan_bprim, scan_frequency, label=folder_labels[i], c=cols[i], lw = 3.0)
        else:
            ax1.plot(scan_bprim, scan_frequency, label=folder_labels[i], lw = 3.0)
        #ax1.set_xlabel("bprim parameter")
        #ax1.set_ylabel("Frequency (normalised GS2 units)")
        ax1.legend(loc='best')

        if use_cols == True:
            ax2.plot(scan_bprim, scan_growth_rate, label=folder_labels[i], c=cols[i], lw = 3.0)
        else:
            ax2.plot(scan_bprim, scan_growth_rate, label=folder_labels[i], lw = 3.0)
        #ax2.set_xlabel("bprim parameter")
        #ax2.set_ylabel("Growth rate (normalised GS2 nuits)")
        ax2.legend(loc='best', prop={'size': legend_font})
    plt.tight_layout()
    plt.show()

    return

def construct_stability_array(stability_df):
    """Construct stability array. Achieve high speed by converting stability_df.stability into a
    (alpha x shat x theta0 array)"""

    theta0 = sorted(set(stability_df.theta0))
    theta0_zero_idx = math.floor(len(theta0)/2)
    print("theta0[theta0_zero_idx] = ", theta0[theta0_zero_idx])

    # sort stability_df based on beta_prime (descending order), then shear (ascending order)
    df2 = stability_df.sort_values(["beta_prime", "shear"], ascending=[False, True])
    #print("df2 = ", df2)

    stability_df_array = (df2.stability.to_numpy())
    stability_df_array.shape = (len(alpha), len(shat), len(theta0))
    #print("stability_df_array = ", stability_df_array)
    #print("")

    for i in range(0, len(alpha)):
        for j in range(0, len(shat)):
            stab_theta0 = stability_df_array[i, j, :]
            #print("theta0, stab_theta0 = ", stab_theta0)

            if stab_theta0[theta0_zero_idx] > 10**-5:
                stability_array[j, i] = 3.0
            else:
                # Check if there are any values of theta0 which are unstable

                #print("i, j = ", i, j)
                if (stab_theta0 > 10**-5).any():
                    stability_array[j, i] = 2.0
                    ## Plot the stability(thtet0) for this simulation
                    #plot_stab_theta0(theta0, stab_theta0, alpha[i], shat[j])
                else:
                    stability_array[j, i] = 1.0

    return [alpha, shat, stability_array]

def plot_omega_for_fprim_tprim_scans(folder_name, normalisations=[], legend_font=8, title=False):
    """Plot a colorplot for growth and frequency, as a function of 2 parameters"""

    # Find the pickled file, load it and extract the values
    omega_dict = load_omega_for_param_scan(folder_name)
    [param1_name, param2_name] = omega_dict["parameter names"]  # Should be fprim, tprim
    param1 = omega_dict[param1_name] ; param2 = omega_dict[param2_name]
    fprim_list = []; tprim_list = []
    for iter in range(len(param1)):
        #print("iter = ", iter)
        #print("param1[iter] = ", param1[iter])
        fprim_val = param1[iter][0][0]
        #print("fprim_val = ", fprim_val)
        fprim_list.append(param1[iter][0][0])
        tprim_list.append(param2[iter][0][0])
    # print("fprim_list = ", fprim_list)
    #print("param1 = ", param1)
    fprim = sorted(set(fprim_list)) # unique values of param1 in ascending order
    tprim = sorted(set(tprim_list)) # unique values of param2 in ascending order
    # print("fprim = ", fprim)
    # print("tprim = ", tprim)
    fprim_array = np.array(fprim); tprim_array = np.array(tprim)
    growth_rate_array = np.zeros((len(fprim), len(tprim)))
    freq_array = np.zeros((len(fprim), len(tprim)))
    #print("param1, param2 = ", param1_set, param2_set)
    print("fprim.shape, tprim.shape, freq_array.shape = ", len(fprim), len(tprim), freq_array.shape)

    for iter in range(len(omega_dict["growth rate"])):
        #print("")
        p1 = param1[iter][0][0]; p2 = param2[iter][0][0]
        idx1 = np.where((abs(fprim_array - p1) < 10**-5))[0]
        idx2 = np.where((abs(tprim_array - p2) < 10**-5))[0]
        #print("idx1, idx2 = ", idx1, idx2)
        growth_rate_array[idx1, idx2] = omega_dict["growth rate"][iter]
        freq_array[idx1, idx2] = omega_dict["frequency"][iter]
        #print("growth_rate = ", omega_dict["growth rate"][iter])

    #print("growth_rate_array = ", growth_rate_array)
    #print("freq_array = ", freq_array)
    #sys.exit()
    fig = plt.figure(figsize=[8, 8])
    ax1 = fig.add_axes([0.15, 0.54, 0.4, 0.4])
    ax2 = fig.add_axes([0.15, 0.1, 0.4, 0.4], sharex=ax1)
    cbax1 = fig.add_axes([0.6, 0.52, 0.03, 0.4])
    cbax2 = fig.add_axes([0.6, 0.05, 0.03, 0.4])
    #ax1.tick_params('x', which='both', labelsize=8, direction="in")
    ax1.grid(True, which="both", axis="both", linestyle="--")
    ax2.grid(True, which="both", axis="both", linestyle="--")

    # Labels
    ax2.set_xlabel(r"fprim"); #ax1.set_xlabel(r"fprim")
    ax2.set_ylabel(r"tprim", rotation=0, labelpad=25)
    ax1.set_ylabel(r"tprim", rotation=0, labelpad=25)


    if title != False:
        fig.suptitle(title, fontsize=16)

    fprim_shift = help_ncdf.shift_axis_for_pcolormesh(fprim)
    tprim_shift = help_ncdf.shift_axis_for_pcolormesh(tprim)

    freqmesh = ax1.pcolormesh(fprim_shift, tprim_shift, freq_array, cmap='Reds')
    growthmesh = ax2.pcolormesh(fprim_shift, tprim_shift, growth_rate_array, cmap='Reds')
    #print("fprim.shape, tprim.shape, freq_array.shape = ", fprim.shape, tprim.shape, freq_array.shape)
    ax1.contour(fprim, tprim, freq_array, 4)#, colors='k')
    ax2.contour(fprim, tprim, growth_rate_array, 4)#, colors='k')
    #fprim, tprim, freq_array
    # cb = plt.colorbar(ax1, cax = cbaxes)
    # plt.colorbar(mappable=colmesh, label="ibs score")
    plt.colorbar(mappable=freqmesh, label=r"$\omega$", cax=cbax1)
    plt.colorbar(mappable=growthmesh, label=r"$\gamma$", cax=cbax2)
    plt.tight_layout()
    plt.show()

    return

def calculate_bprim(input_filetext):
    """ """
    tprim = []
    fprim = []
    dens = []
    temp = []

    #print("filetext = ", filetext)
    # NB these won't work if there are comments etc. in the lines of interest.
    for line in input_filetext.splitlines():
      #print("line = ", line)
      if re.search("beta = ", line):
          print("aha!")
          beta = float(line.split('!')[0].split('=')[1])
          print("beta = ", beta)
      if re.search("tprim", line):
          tprim.append(float(line.split('!')[0].split('=')[1]))
      if re.search("fprim", line):
          fprim.append(float(line.split('!')[0].split('=')[1]))
      if re.search("temp", line):
          temp.append(float(line.split('!')[0].split('=')[1]))
      if re.search("dens", line):
          dens.append(float(line.split('!')[0].split('=')[1]))

    beta_prime = -beta * np.sum(np.array(temp) * np.array(dens) * (np.array(tprim) + np.array(fprim)))
    return beta_prime

def plot_mode_structure_for_param_scan(folder_name):
    """For all simulations in the folder, plot (1) abs(phi) (2) abs(apar)
    and (3) abs(bpar) against theta """

    path = LIN_SIM_PATH + folder_name
    output_file_list = glob.glob(path + '/*.out.nc')
    print("output_file_list = ", output_file_list)

    for nc_filename in output_file_list:
        try:
            [theta, phi] = help_ncdf.extract_data_from_ncdf(nc_filename, "theta", "phi")#, "apar", "bpar")

            phi = phi[0][0]; #apar = apar[0][0]; bpar = bpar[0][0]
            fig1 = plt.figure(figsize=(12, 12))
            ax1 = fig1.add_subplot(111)
            ax1.plot(theta/np.pi, abs(phi)/max(abs(phi)))

            # ax2 = fig1.add_subplot(312, sharex=ax1)
            # ax2.grid(True, which="both", axis="both", linestyle="--")
            # ax2.plot(theta/np.pi, abs(apar)/max(abs(apar)))
            # ax3 = fig1.add_subplot(313, sharex=ax1)
            # ax3.grid(True, which="both", axis="both", linestyle="--")
            # ax3.plot(theta/np.pi, abs(bpar)/max(abs(bpar)))

            ax1.xaxis.set_major_locator(plt.MultipleLocator(2))
            ax1.xaxis.set_minor_locator(plt.MultipleLocator(1))
            ax1.grid(True, which="both", axis="both", linestyle="--")
            ax1.set_yscale("log");# ax2.set_yscale("log"); ax3.set_yscale("log");
            ax1.set_ylabel("abs(phi), normalised"); # ax2.set_ylabel("abs(apar), normalised"); ax3.set_ylabel("abs(bpar), normalised")
            ax1.set_ylim(10**-3, 1.2);# ax2.set_ylim(10**-3, 1.2); ax3.set_ylim(10**-3, 1.2)
            ax1.set_xlabel(r"$\theta/\pi$")

            fig1.suptitle(re.split((folder_name + "/"), nc_filename)[1] +"_mode")
            #plt.show()
            plt.savefig(nc_filename +"_mode.png")
            plt.close()

        except KeyError as err:
            print("nc_filename = ", nc_filename)
            print("KeyError: ", err)

def plot_mode_structure_kxt_for_param_scan(folder_name):
    """For all simulations in the folder, plot (1) abs(phi) (2) abs(apar)
    and (3) abs(bpar) against theta """

    path = LIN_SIM_PATH + folder_name
    output_file_list = glob.glob(path + '/*.out.nc')
    print("output_file_list = ", output_file_list)

    for nc_filename in output_file_list:
        try:
            #print((ncdf2dict.ncdf2dict(nc_filename)).keys())
            [theta, phi, theta0, ky, shat] = help_ncdf.extract_data_from_ncdf(nc_filename, "theta", "phi", "theta0", "ky", "shat")

            kxt = ky*shat*(theta-theta0)[0,:]
            phi = phi[0][0]; #apar = apar[0][0]; bpar = bpar[0][0]
            fig1 = plt.figure(figsize=(12, 12))
            ax1 = fig1.add_subplot(111)
            ax1.grid(True, which="both", axis="both", linestyle="--")
            ax1.plot(kxt, abs(phi)/max(abs(phi)))

            # ax2 = fig1.add_subplot(312, sharex=ax1)
            # ax2.grid(True, which="both", axis="both", linestyle="--")
            # ax2.plot(kxt, abs(apar)/max(abs(apar)))
            # ax3 = fig1.add_subplot(313, sharex=ax1)
            # ax3.grid(True, which="both", axis="both", linestyle="--")
            # ax3.plot(kxt, abs(bpar)/max(abs(bpar)))

            ax1.set_yscale("log"); # ax2.set_yscale("log"); ax3.set_yscale("log");
            ax1.set_ylabel("abs(phi), normalised"); # ax2.set_ylabel("abs(apar), normalised"); ax3.set_ylabel("abs(bpar), normalised")
            ax1.set_ylim(10**-3, 1.2); #ax2.set_ylim(10**-3, 1.2); ax3.set_ylim(10**-3, 1.2)
            ax1.set_xlabel(r"$k_{xt}$")

            fig1.suptitle(re.split((folder_name + "/"), nc_filename)[1] +"_mode")
            #plt.show()
            plt.savefig(nc_filename +"_kxt_mode.png")
            plt.close()

        except KeyError as err:
            print("nc_filename = ", nc_filename)
            print("KeyError: ", err)

def convert_exponent_safe(val_bytes):
    """Convert a bytes string to a float, correctly handling values where the "e"
    is missing from the exponent e.g. 2.55+144"""

    try:
        val_float = float(val_bytes)
    except ValueError:
        #print("Hello. Value error occured, trying to fix")
        #print("val_str = ", val_bytes)
        val_str = val_bytes.decode("utf-8")
        #print("val_str = ", val_str)
        new_val_str = re.sub("\+", "E+", val_str)
        #print("new_val_str = ", new_val_str)
        val_float = float(new_val_str)

    return val_float

def get_abs(re_part, im_part):

    # There's a risk that the real and imaginary
    # components separately are finite, but the abs(value) is not.
    scaling_val = max(re_part)
    re_part = re_part/scaling_val; im_part = im_part/scaling_val
    abs_val = np.sqrt(re_part**2 + im_part**2)
    return abs_val * scaling_val

if __name__ == "__main__":

    print("hello")

    ############################################################################
    ## Open 2 folders, selct the beta = 2% simulation and extract theta
    ############################################################################
    # folder1 = LIN_SIM_PATH + "cmiller_res6/"
    # sim_name1 = folder1 + "cmiller_res6_0.0200"
    # nc_file1 = sim_name1 + ".out.nc"
    # data_dict1 = ncdf2dict.ncdf2dict(nc_file1)
    #
    # folder2 = LIN_SIM_PATH + "cmiller_high_res/"
    # sim_name2 = folder2 + "cmiller_high_res_0.0200"
    # nc_file2 = sim_name2 + ".out.nc"
    # data_dict2 = ncdf2dict.ncdf2dict(nc_file2)
    #
    # print("keys = ", data_dict1.keys())
    # #sys.exit()
    # theta1 = data_dict1["theta"]
    # theta2 = data_dict2["theta"]

    ############################################################################
    ## Code to compare geometrical quantities between the 2 selected simulations
    ############################################################################
    # vars_to_plot = ["bmag", "bpol", "gbdrift", "gbdrift0", "gds2", "gds21", "gds22", "cvdrift", "cvdrift0"]
    # vars_to_plot = ['phi2', 'phi', 'phi_igomega_by_mode', 'phi2_by_mode']
    #                 # "", ""]
    # fig = plt.figure()
    # ax1 = fig.add_subplot(331); ax2 = fig.add_subplot(332); ax3 = fig.add_subplot(333)
    # ax4 = fig.add_subplot(334); ax5 = fig.add_subplot(335); ax6 = fig.add_subplot(336)
    # ax7 = fig.add_subplot(337); ax8 = fig.add_subplot(338); ax9 = fig.add_subplot(339)
    #
    # myaxes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
    #
    # for i, var_name in enumerate(vars_to_plot):
    #     myax = myaxes[i]
    #     var1 = data_dict1[var_name]
    #     var2 = data_dict2[var_name]
    #     if (len(theta1) != len(var1) or len(theta2) != len(var2)):
    #         print("var_name, len(var1), len(theta1) = ", var_name, len(var1), len(theta1))
    #     else:
    #         myax.plot(theta1, var1, label=r"$n_\theta = 256$")
    #         myax.plot(theta2, var2, label=r"$n_\theta = 64$", linestyle="--")
    #         myax.set_xlabel(r"$\theta$")
    #         myax.set_ylabel(var_name)
    #         myax.legend(loc="best")
    #         myax.set_xlim(-16, 16)
    #
    # plt.show()

    ############################################################################
    ## Code to examine mode structures (phi(theta)) of simulations
    ############################################################################

    # phi_sim1 = data_dict1["phi"][0][0]
    # phi_sim2 = data_dict2["phi"][0][0]
    #
    # #### Normalise the mode structures, so the maximum amplitude is 1
    # print("max(phi_sim1) = ", max(phi_sim1))
    # phi_sim1 /= max(abs(phi_sim1))
    # print("max(phi_sim1) = ", max(phi_sim1))
    # phi_sim2 /= max(abs(phi_sim2))
    # plt.figure()
    # plt.plot(theta1, abs(phi_sim1), label=r"$n_\theta = 256$")
    # plt.plot(theta2, abs(phi_sim2), label=r"$n_\theta = 64$", linestyle="--")
    # plt.yscale("log")
    # plt.legend(loc="best")
    # plt.show()

    ############################################################################
    ## Code to plot and save mode structure for many different simulations
    ############################################################################

    # plot_upar_function_vs_theta_ntheta("ntot", "ntot_ions_", spec=0)
    # plot_upar_function_vs_theta_ntheta("ntot", "ntot_electrons_", spec=1)
    # plot_upar_function_vs_theta_ntheta("density", "dens_ions_", spec=0)
    # plot_upar_function_vs_theta_ntheta("density", "dens_electrons_", spec=1)
    # plot_upar_function_vs_theta_ntheta("upar", "upar_ions_", spec=0)
    # plot_upar_function_vs_theta_ntheta("upar", "upar_electrons_", spec=1)
    # plot_upar_function_vs_theta_ntheta("apar", "apar_")
    # plot_upar_function_vs_theta_ntheta("bpar", "bpar_")
    # plot_upar_function_vs_theta_ntheta("phi", "phi_")

    ############################################################################
    ## Code to scatter growth rate vs frequency for different simulations
    ############################################################################

    # (param_list1, frequency_list1, growth_rate_list1,
    #         freq_errors1, growth_rate_errors1) = calculate_omega_for_ntheta_vals("cmiller_beta3")
    #
    # (param_list2, frequency_list2, growth_rate_list2,
    #         freq_errors2, growth_rate_errors2) = calculate_omega_for_ntheta_vals_4("cmiller_beta3_opt_source")
    #
    # marker_size = 50.
    # label_info = r"opt_source"; mytitle="comparison with opt_sorce=true"
    # plt.figure()
    # plt.scatter(frequency_list1[0], growth_rate_list1[0], c="green", s=marker_size, marker="x", label=r"$n_\theta=64$")
    # plt.scatter(frequency_list1[1], growth_rate_list1[1], c="orange", s=marker_size, marker="x", label=r"$n_\theta=128$")
    # plt.scatter(frequency_list1[2], growth_rate_list1[2], c="blue", s=marker_size, marker="x", label=r"$n_\theta=256$")
    # plt.errorbar(frequency_list1, growth_rate_list1, xerr=freq_errors1, yerr=growth_rate_errors1, fmt='none',
    # c='r', capsize=3)
    # plt.scatter(frequency_list2[0], growth_rate_list2[0], c="black", marker="o", s=marker_size, label=r"$n_\theta=32$, {0}".format(label_info))
    # plt.scatter(frequency_list2[1], growth_rate_list2[1], c="green", marker="o", s=marker_size, label=r"$n_\theta=64$, {0}".format(label_info))
    # plt.scatter(frequency_list2[2], growth_rate_list2[2], c="orange", marker="o", s=marker_size, label=r"$n_\theta=128$, {0}".format(label_info))
    # plt.scatter(frequency_list2[3], growth_rate_list2[3], c="blue", marker="o", s=marker_size, label=r"$n_\theta=256$, {0}".format(label_info))
    # plt.errorbar(frequency_list2, growth_rate_list2, xerr=freq_errors2, yerr=growth_rate_errors2, fmt='none',
    # c='r', capsize=3)
    # plt.legend(loc="best")
    # plt.xlabel(r"$\omega$")
    # plt.ylabel(r"$\gamma$")
    # plt.title(mytitle)
    # plt.show()



    ############################################################################

    #make_folders_and_process_scan("gcsp_highres_ky_beta_lowky")

    #help_ncdf.view_ncdf_variables('cilr_aky_010.out.nc')
    # #compare_simulations_different_machines()
    #save_nonlinear_phi2()



    # compare_simulations(["nimpurity_m7q7_diags_0.0100.out.nc", "nimpurity_m7q7_diags_0.0100_nperiod5.out.nc"],
    #                      title="m=7, q=7, beta=1%",
    #                     labels=["nperiod=3", "nperiod=5"])
