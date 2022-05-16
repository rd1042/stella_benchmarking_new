""" """
import sys
sys.path.append("../postprocessing_tools")
from omega_calculator import write_gs2_data_to_plaintext
import glob
import re

gs2_fapar0_fbpar1_me1_folder = "gs2_fapar0_fbpar1_me1_beta_scan"
gs2_fapar1_fbpar0_me1_folder = "gs2_fapar1_fbpar0_me1_beta_scan"
gs2_fapar1_fbpar1_me1_folder = "gs2_fapar1_fbpar1_me1_beta_scan"

def save_gs2_sims_in_folder(folder_longname):
    """Find all GS2 data in a folder and save it"""
    outnc_file_list = glob.glob(folder_longname + "/*.out.nc")

    for outnc_longname in outnc_file_list:
        sim_longname = re.split(".out.nc", outnc_longname)[0]
        write_gs2_data_to_plaintext(sim_longname)


if __name__ == "__main__":
    print("hello world")
    #folder_longname = "gs2_fapar1_fbpar0_me1_beta_scan"
    # sim_longname = "gs2_fapar1_fbpar0_me1_beta_scan/beta_0.00800"
    # write_gs2_data_to_plaintext(sim_longname)
    save_gs2_sims_in_folder(gs2_fapar0_fbpar1_me1_folder)
    save_gs2_sims_in_folder(gs2_fapar1_fbpar0_me1_folder)
    save_gs2_sims_in_folder(gs2_fapar1_fbpar1_me1_folder)
