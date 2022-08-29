"""Make simulations which scan vpa, mu, ky. To test field solve
in h and gbar """

import numpy as np
import re
import sys
import os
import shutil
import sys
sys.path.append("../generate_sims/")
import make_param_scans as mps

phi_bpar_h_vpa_scan_folder = "stella_phi_bpar_h_nvpa_vpamax_scan"
phi_bpar_h_vperp_scan_folder = "stella_phi_bpar_h_nmu_vperpmax_scan"
phi_bpar_gbar_vpa_scan_folder = "stella_phi_bpar_gbar_nvpa_vpamax_scan"
phi_bpar_gbar_vperp_scan_folder = "stella_phi_bpar_gbar_nmu_vperpmax_scan"

def construct_nvpa_vpamax_scan(folder_shortname):
    folder_name = "sims/" + folder_shortname
    template_name = folder_name + "/template_input_file"
    template_runscript = None

    vpamax_vals = np.array([0.1, 0.3, 0.5, 0.8, 1, 1.5, 2, 3])
    nvpa_vals = np.array([4, 6, 8, 10, 12, 15, 18, 21, 24, 30, 36, 48, 60, 72, 96, 120, 144])

    # Check folder exists
    if not os.path.isdir(folder_name):
        sys.exit("folder " + folder_name  + " does not exist!")

    # Check if the template file exists
    if not os.path.isfile(template_name):
        sys.exit("template file does not exist!")

    mps.make_input_files_2d(template_name, template_runscript, "vpa_max", "nvgrid",
                vpamax_vals, nvpa_vals, "vpamax_nvpa", folder=folder_name, dtypes=["float", "int"])
    print("folder contents: " + folder_name)
    os.system('ls -ltr ' + folder_name)
    return

def construct_nmu_vperpmax_scan(folder_shortname):
    folder_name = "sims/" + folder_shortname
    template_name = folder_name + "/template_input_file"
    template_runscript = None

    vperpmax_vals = np.array([0.1, 0.3, 0.5, 0.8, 1, 1.5, 2, 3])
    nmu_vals = np.array([4, 6, 8, 10, 12, 15, 18, 21, 24, 30, 36, 48, 60, 72, 96, 120, 144])

    # Check folder exists
    if not os.path.isdir(folder_name):
        sys.exit("folder " + folder_name  + " does not exist!")

    # Check if the template file exists
    if not os.path.isfile(template_name):
        sys.exit("template file does not exist!")

    mps.make_input_files_2d(template_name, template_runscript, "vperp_max", "nmu",
                vperpmax_vals, nmu_vals, "vperpmax_nmu", folder=folder_name, dtypes=["float", "int"])
    print("folder contents: " + folder_name)
    os.system('ls -ltr ' + folder_name)
    return

def make_input_file(template, template_runscript, parameters, values, beware_integers=False, filename_prefix="",
                    filename=False, folder="."):
    """Constructs an input file based on the template, changing an arbitrary number of
    parameters."""

    ## Check if same number of parameters as values
    if len(parameters) != len(values):
        print("Error! len(parameters), len(values) = ", len(parameters), len(values))
        sys.exit()
    template_file = open(template, 'r+')
    filetext = template_file.read()
    template_file.close()
    values_str = ""
    for i, parameter in enumerate(parameters):
      value = values[i]
      if beware_integers:
          if isinstance(value, int):
              value_str = str(value)
          elif isinstance(value, float):
            value_str = ('%.4f' % value)
          else:
              print("Error! Type not recognised. value = ", value)
      else:
          value_str = ('%.5f' % value)

      values_str += ("_" + value_str)

      if type(parameter) is str:
          search_string = parameter + " = .*"
          replace_string = parameter + " = " + value_str
          filetext = re.sub(search_string, replace_string, filetext)
      elif type(parameter) is list:
          for each_parameter in parameter:
              search_string = each_parameter + " = .*"
              replace_string = each_parameter + " = " + value_str
              filetext = re.sub(search_string, replace_string, filetext)
      else:
          print("type(parameter) not recognised! Aborting" )
          print("type(parameter) = ", type(parameter))
          sys.exit()

    if filename != False:
        new_filename = filename
    else:
        #print("folder = ", folder)
        infile_shortname = filename_prefix + values_str + ".in"
        new_filename = folder +  "/" + infile_shortname

    #print("new_filename = ", new_filename)
    #sys.exit()
    new_file = open(new_filename, "w+")
    new_file.write(filetext)
    new_file.close()

    ## Open and modify run script.
    template_file = open(template_runscript, 'r+')
    filetext = template_file.read()
    template_file.close()
    filetext = re.sub("input.in", infile_shortname, filetext)
    new_runscript_name = folder + "/run_beta" + values_str + ".sh"
    new_file = open(new_runscript_name, "w+")
    new_file.write(filetext)
    new_file.close()
    return

def make_sims_phi_bpar_test():
    """ """
    #construct_nvpa_vpamax_scan(phi_bpar_h_vpa_scan_folder)
    construct_nmu_vperpmax_scan(phi_bpar_h_vperp_scan_folder)

    return

if __name__ == "__main__":
    make_sims_phi_bpar_test()
