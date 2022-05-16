""" """
import sys
sys.path.append("../postprocessing_tools")
from omega_calculator import write_gs2_data_to_plaintext
import glob
import re

# gs2
gs2_folder = "gs2_kykx_scan"
gs2_me1_folder = "gs2_me1_kykx_scan"
gs2_adiabatic_folder = "gs2_adiabatic_kykx_scan"

def save_gs2_sims_in_folder(folder_longname):
    """Find all GS2 data in a folder and save it"""
    outnc_file_list = glob.glob(folder_longname + "/*.out.nc")

    for outnc_longname in outnc_file_list:
        sim_longname = re.split(".out.nc", outnc_longname)[0]
        write_gs2_data_to_plaintext(sim_longname)


if __name__ == "__main__":
    print("hello world")
    save_gs2_sims_in_folder(gs2_folder)
    save_gs2_sims_in_folder(gs2_me1_folder)
    save_gs2_sims_in_folder(gs2_adiabatic_folder)
