""" """

import numpy as np
import pickle
import matplotlib.pyplot as plt
import math
import sys

from extract_sim_data import get_omega_data, get_phiz_data, get_aparz_data, get_bparz_data


folder1 = "sims/w003_es_linear_master_src_h_comparison/"
master_explicit_outnc = folder1 + "master_explicit.out.nc"
master_str_m_implicit_outnc = folder1 + "master_str_m_implicit.out.nc"
src_h_explicit_outnc = folder1 + "src_h_explicit.out.nc"
src_h_str_m_implicit_outnc = folder1 + "src_h_str_m_implicit.out.nc"

def make_comparison_plots_for_thesis():
    """ """
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    outnc_longnames = [master_explicit_outnc,
                    master_str_m_implicit_outnc,
                    src_h_explicit_outnc,
                    src_h_str_m_implicit_outnc]
    for counter, outnc_longname in enumerate(outnc_longnames):



    return

if __name__ == "__main__":
    print("Hello world")
