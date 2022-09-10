""" """

import sys
sys.path.append("../postprocessing_tools")
import matplotlib.pyplot as plt
import numpy as np
import glob
import re
import pickle

IMAGE_DIR = "./images/"

leapfrog_nonlinear_outnc = "stella_nonlinear_adiabatic_leapfrog/input.out.nc"
master_nonlinear_outnc = "stella_nonlinear_adiabatic_master/input.out.nc"


def plot_phi2_t(pickle_longnames):
    """ """

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    for pickle_longname in pickle_longnames:
        shortname = re.split("/", pickle_longname)[-2]
        myfile = open(pickle_longname, "rb")
        [t, kx, ky, phi2, phi2_kxky] = pickle.load(myfile)
        myfile.close()
        ax1.plot(t, phi2, label=shortname, lw=2)

    ax1.set_yscale("log")
    ax1.set_xlabel(r"$t$")
    ax1.set_ylabel(r"phi2")
    ax1.legend(loc="best")
    plt.show()

    return

def plot_phi2_t_all():
    """ """
    pickle_longnames = glob.glob("sims/*/*.summary_pickle")
    plot_phi2_t(pickle_longnames)

    return

if __name__ == "__main__":
    print("Hello world")
    plot_phi2_t_all()
