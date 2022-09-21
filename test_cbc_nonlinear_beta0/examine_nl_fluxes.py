""" """

import sys
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pickle
import re
from scipy.interpolate import griddata

def extract_data_from_ncdf_with_xarray(sim_name, *args):
    """Extract data arrays from the NetCDF file for a simulation. Extracts all
    data given in *args. TO IMPLEMENT:
     - Tell you if the file doens't exist
     - If args don't exist, view_ncdf variables"""

    # Convert the input file to a python dictionary
    data = xr.open_dataset(sim_name)
    datalist = []

    # Extract each array specified in args, and add it to the list
    for arg in args:
        datalist.append(data[arg])

    return datalist

if __name__ == "__main__":
    outnc_longname = "sims/test_stella_master_with_fluxes/input.out.nc"
    [t, kx, ky, z, phi2_t, phi2_tkxky, pflx_kxky, vflx_kxky, qflx_kxky] = extract_data_from_ncdf_with_xarray(outnc_longname,
                                    "t", "kx", "ky", "zed", "phi2", "phi2_vs_kxky",
                                    "pflx_kxky", "vflx_kxky", "qflx_kxky")
    print("qflx_kxky = ", qflx_kxky)
    t = np.array(t)
    z = np.array(z)
    kx = np.array(kx)
    ky = np.array(ky)
    phi2_t = np.array(phi2_t)
    phi2_tkxky = np.array(phi2_tkxky)
    qflx_kxky = np.array(qflx_kxky)

    print("len(t) = ", len(t))
    print("len(kx) = ", len(kx))
    print("len(ky) = ", len(ky))
    print("len(z) = ", len(z))
    print("phi2_tkxky.shape = ", phi2_tkxky.shape)
    print("qflx_kxky.shape = ", qflx_kxky.shape) # t, spec, tube, z, kx, ky

    t_idxs = [0, int(0.25*len(t)), int(0.5*len(t)), int(0.75*len(t)), -1]
    z_idx = int(len(z)/2)
    print("t_idxs = ", t_idxs)
    parsed_qflx = qflx_kxky[t_idxs,:,0,z_idx,:,:]
    print("parsed_qflx = ", parsed_qflx)
    print("parsed_qflx.shape = ", parsed_qflx.shape)
