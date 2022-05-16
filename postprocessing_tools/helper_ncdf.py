"""Helper module, containing often-used and mature methods. The aim is for this
script to be thoroughly tested, and then called by other scripts as needed.

Author: Bob Davies"""

import ncdf2dict
import sys
import xarray as xr


def extract_data_from_ncdf(sim_name, *args):
    """Extract data arrays from the NetCDF file for a simulation. Extracts all
    data given in *args. TO IMPLEMENT:
     - Tell you if the file doens't exist
     - If args don't exist, view_ncdf variables"""


    # Convert the input file to a python dictionary
    data = ncdf2dict.ncdf2dict(sim_name)
    datalist = []
    # Extract each array specified in args, and add it to the list
    for arg in args:
        datalist.append(data[arg])

    return datalist

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

def view_ncdf_variables(sim_name):
    """View the names all variables in the netcdf file"""
    data = ncdf2dict.ncdf2dict(sim_name)
    print(list(data.keys()))
    return
