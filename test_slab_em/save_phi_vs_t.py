""" """

import sys
sys.path.append("../postprocessing_tools")
import glob
import re
import numpy as np
from helper_ncdf import extract_data_from_ncdf
import pickle

def save_sim_phi_vs_t(sim_longname):
    """ """
    outnc_longname = sim_longname + ".out.nc"
    t, z, phi_vs_t = extract_data_from_ncdf(outnc_longname, "t", 'zed', 'phi_vs_t')
    # print("phi_vs_t.shape = ", phi_vs_t.shape)  # [n_time, 1 , n_z, 1, 1]
    phi_vs_t = phi_vs_t[:,0,:,0,0]
    outfile_name = sim_longname + ".pickle"
    outfile = open(outfile_name, "wb")
    pickle.dump([z, t, phi_vs_t[-1,:], phi_vs_t[:,int(len(z)*0.5)]], outfile)
    # outfile.write("z \n")
    # outfile.write(str(list(z)))
    # outfile.write("\n t \n")
    # outfile.write(str(list(t)))
    # outfile.write("\n phi_vs_t[-1,:] \n")
    # outfile.write(str(list(phi_vs_t[-1,:])))
    # outfile.write("\n phi_vs_t[:,int(len(z)*0.5)] \n")
    # outfile.write(str(list(phi_vs_t[:,int(len(z)*0.5)])))
    outfile.close()

    return

def save_phi_vs_t_mandell_sims():
    """ """
    infile_longnames = glob.glob("mandell_sims/*.out.nc")
    print("infile_longnames = ", infile_longnames)
    for infile_longname in infile_longnames:
        sim_longname = re.split(".out.nc", infile_longname)[0]
        print("sim_longname = ", sim_longname)
        save_sim_phi_vs_t(sim_longname)
    return

def save_phi_vs_t_mandell_sims2():
    """ """
    infile_longnames = glob.glob("mandell_sims_supermassive_ions/*.out.nc")
    print("infile_longnames = ", infile_longnames)
    for infile_longname in infile_longnames:
        sim_longname = re.split(".out.nc", infile_longname)[0]
        print("sim_longname = ", sim_longname)
        save_sim_phi_vs_t(sim_longname)
    return



if __name__ == "__main__":
    #save_phi_vs_t_mandell_sims()
    save_phi_vs_t_mandell_sims2()
