""" """
import sys
import re
import numpy as np
from numpy import linalg

sim_folder = "sims/stella_ky0.05_stability_tests/"

src_h_fphi1_fapar1_fbpar1_ky005_longname = sim_folder + "input_implicit_ky0.05_src_h_fphi1_fapar1_fbpar1"
src_h_fphi1_fapar1_fbpar1_ky05_longname = sim_folder + "input_implicit_ky0.5_src_h_fphi1_fapar1_fbpar1"
src_h_fphi1_fapar1_fbpar1_ky05_dt01_longname = sim_folder + "input_implicit_ky0.5_src_h_fphi1_fapar1_fbpar1_dt0.1"
src_h_fphi1_fapar0_fbpar0_ky005_longname = sim_folder + "input_implicit_ky0.05_src_h_fphi1_fapar0_fbpar0"
src_h_fphi0_fapar1_fbpar0_ky005_longname = sim_folder + "input_implicit_ky0.05_src_h_fphi0_fapar1_fbpar0"
src_h_fphi0_fapar1_fbpar0_ky05_longname = sim_folder + "input_implicit_ky0.5_src_h_fphi0_fapar1_fbpar0"

def get_response_matrix(sim_filename):
    """ """
    response_longame = sim_filename + ".streaming_response"
    myfile = open(response_longame, "r")
    filetext = myfile.read().strip()
    myfile.close()
    matrices = re.split("matrix size = ", filetext)[1:]
    n_matrices = len(matrices)
    if n_matrices != 1:
        print("Error! n_matrices != 1")
        print("n_matrices = ", n_matrices)
        sys.exit()
    lines = re.split("\n", matrices[0])
    size_str = lines[0].strip()
    size_numbers = re.split("\s+", size_str)
    matrix_dim = int(size_numbers[0])
    print("matrix_dim = ", matrix_dim)
    response_matrix = np.zeros((matrix_dim, matrix_dim))

    for row_idx, line in enumerate(lines[1:]):
        matrix_str = line.strip()
        matrix_entries_str = re.split("\s+", matrix_str)
        for col_idx, entry_str in enumerate(matrix_entries_str):
            entry_str = entry_str[1:-1]
            real_str, im_str = re.split(",", entry_str)
            real_entry = float(real_str)
            response_matrix[row_idx, col_idx] = real_entry

    return response_matrix

def get_condition_number(response_matrix):
    """ """
    condition_number = linalg.cond(response_matrix)
    print("condition_number = ", condition_number)
    # print("norm = ", linalg.norm(response_matrix))
    inversion = linalg.inv(response_matrix)
    # print("norm(inv) = ", linalg.norm(inversion))
    return condition_number

def process_sim(sim_longname):
    """ """
    response_matrix = get_response_matrix(sim_longname)
    print("got response matrix")
     # get_condition_number(response_matrix[:40,:40]) returns inf.
     # get_condition_number(response_matrix) times out.
    return get_condition_number(response_matrix)

if __name__ == "__main__":
    print("Hello world")
    condition_numbers = []
    condition_numbers.append(process_sim(src_h_fphi1_fapar0_fbpar0_ky005_longname))
    condition_numbers.append(process_sim(src_h_fphi0_fapar1_fbpar0_ky005_longname))
    condition_numbers.append(process_sim(src_h_fphi0_fapar1_fbpar0_ky05_longname))
    condition_numbers.append(process_sim(src_h_fphi1_fapar1_fbpar1_ky005_longname))
    condition_numbers.append(process_sim(src_h_fphi1_fapar1_fbpar1_ky05_longname))
    condition_numbers.append(process_sim(src_h_fphi1_fapar1_fbpar1_ky05_dt01_longname))
    print("condition_numbers = ", condition_numbers)
