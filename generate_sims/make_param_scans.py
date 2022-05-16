"""Code to generate scans in beta' for both stella (running on MARCONI)
and GS2 (running on VIKING) """

import numpy as np
import re
import sys
import os
import shutil

def construct_beta_scan(folder_name, value_list, sim_type):
    template_name = folder_name + "/template_input_file"
    template_runscript = folder_name + "/template_run_script.sh"

    # Check folder exists
    if not os.path.isdir(folder_name):
        sys.exit("folder " + folder_name  + " does not exist!")

    # Check if the template file exists
    if not os.path.isfile(template_name):
        sys.exit("template file does not exist!")


    make_input_files_1d(template_name, template_runscript, 'beta', value_list, 'beta', folder=folder_name)
    shutil.copy('../templates_and_scripts/submit_jobs.sh', folder_name); os.chmod((folder_name + '/submit_jobs.sh'), 0o777)
    print("folder contents: " + folder_name)
    os.system('ls -ltr ' + folder_name)
    return

def construct_ky_kx_scan(folder_name, ky_vals, kx_vals, sim_type):
    template_name = folder_name + "/template_input_file"
    template_runscript = folder_name + "/template_run_script.sh"

    # Check folder exists
    if not os.path.isdir(folder_name):
        sys.exit("folder " + folder_name  + " does not exist!")

    # Check if the template file exists
    if not os.path.isfile(template_name):
        sys.exit("template file does not exist!")

    if sim_type == "stella":
        construct_ky_kx_scan_stella(template_name, template_runscript, ky_vals, kx_vals, folder=folder_name)
    elif sim_type == "gs2":
        make_input_files_2d(template_name, template_runscript, "aky", "akx", ky_vals, kx_vals, "ky_kx", folder=folder_name)
    shutil.copy('../templates_and_scripts/submit_jobs.sh', folder_name); os.chmod((folder_name + '/submit_jobs.sh'), 0o777)
    print("folder contents: " + folder_name)
    os.system('ls -ltr ' + folder_name)
    return

def construct_ky_kx_scan_stella(template_name, template_runscript, ky_vals, kx_vals, folder):
    """ """

    for ky_val in ky_vals:
        for kx_val in kx_vals:
            make_input_file(template_name, template_runscript,
                    ["aky_min", "aky_max", "akx_min", "akx_max"],
                    [ky_val, ky_val, kx_val, kx_val], filename_prefix="ky_kx",
                        folder=folder)

    return

def make_input_files_1d(template, template_runscript, parameter, value_list, filename_prefix, folder="."):
  """Constructs many input files, changing a single parameter to values based
  a list"""

  for value in value_list:
      make_input_file(template, template_runscript, [parameter], [value], filename_prefix=filename_prefix,
                        folder=folder)

  return

def make_input_files_2d(template, template_runscript, parameter1, parameter2, value_list1,
                        value_list2, filename_prefix, folder="."):
  """Constructs many input files, changing a single parameter to values based
  a list"""

  for value1 in value_list1:
      for value2 in value_list2:
          make_input_file(template, template_runscript, [parameter1, parameter2], [value1, value2], filename_prefix=filename_prefix,
                        folder=folder)

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
          value_str = ('%.8f' % value)

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

    if template_runscript is not None:
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

if __name__ == "__main__":
    print("Hello world")
