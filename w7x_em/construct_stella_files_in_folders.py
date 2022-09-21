"""Script to create multiple input files based on a template, in order to
perform parameter scans in GS2"""

import re
import numpy as np
import sys
import os
import shutil

def make_input_file(template, folder_longname, parameters, values, filename=None):
    """Constructs an input file based on the template, changing an arbitrary number of
    parameters."""

    template_file = open(template, 'r+')
    filetext = template_file.read()
    template_file.close()
    values_str = ""
    for i, parameter in enumerate(parameters):
      value = values[i]
      value_str = ('%.4f' % value)
      values_str += ("_" + value_str)
      #print("values_str = ", values_str)
      # Make a copy of the template with new name
      #shutil.copyfile(template, new_filename)

      # Open the new file
      ### Sometimes this is "<parameter> = <value>", sometimes "<parameter>=<value>"
      ### Distinguish.
      search_string = parameter + "=.*"
      replace_string = parameter + "=" + value_str
      #print("search_string = ", search_string)
      #print("replace_string = ", replace_string)
      filetext = re.sub(search_string, replace_string, filetext)
      #print("filetext = ", filetext)

    if filename is not None:
      new_file_shortname = filename
    else:
      new_file_shortname = "input.in"
    print("folder_longname, new_file_shortname = ", folder_longname, new_file_shortname)
    file_longname = folder_longname + "/" + new_file_shortname
    new_file = open(file_longname, "w+")
    new_file.write(filetext)
    new_file.close()


    return

def make_ky_scan(parent_folder_name, ky_vals):
    """ """
    ## Look for a template_input_file and a template_runscript
    template = parent_folder_name + "/template_input_file"
    template_runscript = parent_folder_name + "/template_runscript"
    if not os.path.exists(template):
        print("template doesn't exist! Exitting")
        sys.exit()
    if not os.path.exists(template_runscript):
        print("template_runscript doesn't exist! Exitting")
        sys.exit()

    ## for each ky val, make an appropriately named folder and add the input file
    ## and runscript
    for ky_val in ky_vals:
        folder_shortname = "ky_{:.4f}".format(ky_val)
        folder_longname = parent_folder_name + "/" + folder_shortname
        os.mkdir(folder_longname)
        runscript_longname = folder_longname + "/run_stella.sh"
        make_input_file(template, folder_longname, ["aky_min", "aky_max"], [ky_val, ky_val])

        shutil.copy(template_runscript, runscript_longname)

    template_runscript

    shutil.copy('../templates_and_scripts/submit_jobs.sh', parent_folder_name); os.chmod((parent_folder_name + '/submit_jobs.sh'), 0o777)

    return

def make_fprim_tprim_ky_scan(parent_folder_name, fprim_vals, tprim_vals, ky_vals):
    """ """
    ## Look for a template_input_file and a template_runscript
    template = parent_folder_name + "/template_input_file"
    template_runscript = parent_folder_name + "/template_runscript"
    if not os.path.exists(template):
        print("template doesn't exist! Exitting")
        sys.exit()
    if not os.path.exists(template_runscript):
        print("template_runscript doesn't exist! Exitting")
        sys.exit()

    ## for each ky val, make an appropriately named folder and add the input file
    ## and runscript
    for fprim_val in fprim_vals:
        for tprim_val in tprim_vals:
            for ky_val in ky_vals:
                folder_shortname = "fprim_{:.4f}_tprim_{:.4f}_ky_{:.4f}".format(fprim_val, tprim_val, ky_val)
                folder_longname = parent_folder_name + "/" + folder_shortname
                os.mkdir(folder_longname)
                make_input_file(template, folder_longname, ["fprim", "tprim", "aky_min", "aky_max"],
                                [fprim_val, tprim_val, ky_val, ky_val])
                runscript_longname = folder_longname + "/run_stella.sh"
                shutil.copy(template_runscript, runscript_longname)

    template_runscript

    shutil.copy('../templates_and_scripts/submit_jobs.sh', parent_folder_name); os.chmod((parent_folder_name + '/submit_jobs.sh'), 0o777)

    return

def make_beta_scan(parent_folder_name, beta_vals):
    """ """
    ## Look for a template_input_file and a template_runscript
    template = parent_folder_name + "/template_input_file"
    template_runscript = parent_folder_name + "/template_runscript"
    if not os.path.exists(template):
        print("template doesn't exist! Exitting")
        sys.exit()
    if not os.path.exists(template_runscript):
        print("template_runscript doesn't exist! Exitting")
        sys.exit()

    ## for each ky val, make an appropriately named folder and add the input file
    ## and runscript
    for beta_val in beta_vals:
        folder_shortname = "beta_{:.4f}".format(beta_val)
        folder_longname = parent_folder_name + "/" + folder_shortname
        os.mkdir(folder_longname)
        runscript_longname = folder_longname + "/run_stella.sh"
        make_input_file(template, folder_longname, ["beta"], [beta_val])

        shutil.copy(template_runscript, runscript_longname)

    shutil.copy('../templates_and_scripts/submit_jobs.sh', parent_folder_name); os.chmod((parent_folder_name + '/submit_jobs.sh'), 0o777)

    return

# if __name__ == "__main__":
#     print("Hello world")
#     fprim_vals = np.arange(0, 8.5, 0.5)
#     tprim_vals = np.arange(0, 8.5, 0.5)
#     folder_name = "../fprim_tprim_scan_theta0_piby4"
#     make_input_files_2d(folder_name, ["fprim", "tprim"], [fprim_vals, tprim_vals], (folder_name + "/fprim_tprim"))
