""" """

import re
import numpy as np

fapar1_fbpar0_centered_dgdz_output_name = "sims/test_write_size_of_streaming_terms/fapar1_results.txt"
fapar1_fbpar0_output_name = "sims/test_write_size_of_streaming_terms/for_ref_fapar1_results.txt"
fapar0_fbpar0_output_name = "sims/test_write_size_of_streaming_terms/for_ref_fapar0_results.txt"

def plot_terms(plaintext_longname):
    """ """
    myfile = open(plaintext_longname, "r")
    filetext = myfile.read()
    vpa_vals = []
    mu_vals = []
    streaming_term_vals = []

    for line in re.split("\n", filetext.strip()):
        stripped_line = line.strip()
        [text, vals] = re.split("=", stripped_line)
        [vpa_str, mu_str, stream_term_str] = re.split("\s+", vals.strip())
        vpa_vals.append(float(vpa_str))
        mu_vals.append(float(mu_str))
        streaming_term_vals.append(float(stream_term_str))

    vpa_array = np.array(vpa_vals)
    mu_array = np.array(mu_vals)
    streaming_term_array = np.array(streaming_term_vals)

    print("np.mean(streaming_term_array) = ", np.mean(streaming_term_array))

    return

def compare_terms_for_sims():
    """ """
    plot_terms(fapar0_fbpar0_output_name)
    plot_terms(fapar1_fbpar0_output_name)
    plot_terms(fapar1_fbpar0_centered_dgdz_output_name)

    return

if __name__ == "__main__":
    compare_terms_for_sims()
