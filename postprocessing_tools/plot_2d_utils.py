"""Some useful functions for making contour plots, colorplots etc. of data z(x, y) Much of the
difficulty in this is around the form of z, x and y (meshgrids, unique values, 1D arrays etc.)
To this end we define 3 different types of data which are commonly used:
 - meshgrid: a regular set of coordinates in x or y (regular 2D array)
 - unique_list: for regularly grid of x and y coordinates, the unique values of x and y
 - long_list: a list of x, y or z values, with as many entries as there are data points.
 Useful if the data is not regular in x, y. NB long_lists of x, y data is in general
 non-unique and non-regular."""

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

def unique_arrays2longlists(unique_colvals, unique_rowvals, data_arrays):
    """The rowvals are turned into the y coordinates, colvals into x coordinates.
    A little confusing because sometimes, the data array is given as [x,y], so in the array
    x is the row val."""

    xval_longlist = []; yval_longlist = []; data_longlists = []
    #print("len(unique_colvals), len(unique_rowvals), data_arrays[0].shape = ", len(unique_colvals), len(unique_rowvals), data_arrays[0].shape )


    for yval in unique_rowvals:
        for xval in unique_colvals:
            xval_longlist.append(xval)
            yval_longlist.append(yval)

    for data_array in data_arrays:

        data_longlist = []
        if len(unique_colvals) == len(data_array[:,0]):
            for row_idx in range(0, len(unique_rowvals)):
                for col_idx in range(0, len(unique_colvals)):
                    data_longlist.append(data_array[col_idx, row_idx])

            data_longlists.append(data_longlist)
        else:
            print("Trying the other way around")
            for row_idx in range(0, len(unique_rowvals)):
                for col_idx in range(0, len(unique_colvals)):
                    data_longlist.append(data_array[row_idx, col_idx])

            data_longlists.append(data_longlist)

    return xval_longlist, yval_longlist, data_longlists

def construct_meshgrid(xval_longlist, yval_longlist, data_longlists, n_xpoints=100, n_ypoints=100,
                       interp_kind="nearest", fill_val="NaN"):
    """Convert x and y longlists into meshgrids, and interpolate data longlists
    onto the meshgrids, with user-defined interpolation method. """

    [y_min, y_max] = [min(yval_longlist), max(yval_longlist)]
    [x_min, x_max] = [min(xval_longlist), max(xval_longlist)]

    # Define new fprim, tprim coordinates
    y_mesh = np.linspace(y_min, y_max, n_ypoints)
    x_mesh = np.linspace(x_min, x_max, n_ypoints)
    #print("x_mesh = ", x_mesh)

    # Turn into a "meshgrid"
    x_mesh, y_mesh = np.meshgrid(x_mesh, y_mesh)

    interpolated_data = []  # To store output.

    for data_longlist in data_longlists:
        # Interpolate values onto new grid.
        # print("xval_longlist = ", xval_longlist)
        # print("yval_longlist = ", yval_longlist)
        # print("data_longlist = ", data_longlist)
        # #print("data_longlist.shape() = ", data_longlist.shape)
        # print("len(xval_longlist), len(yval_longlist) = ", len(xval_longlist), len(yval_longlist) )
        interpolated_data.append(interp_single_onto_meshgrid(xval_longlist, yval_longlist,
            np.asarray(data_longlist), x_mesh, y_mesh, interp_kind=interp_kind, fill_val=fill_val))
        # interpolated_data.append(griddata((xval_longlist, yval_longlist),
        #                   np.asarray(data_longlist), (x_mesh, y_mesh),
        #                   method=interp_kind, fill_value=100.0))

    return x_mesh, y_mesh, interpolated_data


def interp_many_onto_meshgrid(xval_longlist, yval_longlist, data_longlists, x_meshrgid, y_meshgrid,
                         interp_kind="nearest", fill_val="NaN"):

    data_meshgrids = []
    for data_longlist in data_longlists:
        #data_longlist = np.ones(len(xval_longlist))
        data_meshgrids.append(griddata((xval_longlist, yval_longlist),np.asarray(data_longlist),(x_meshrgid, y_meshgrid),
                         method=interp_kind, fill_value=fill_val))
    return data_meshgrids

def interp_single_onto_meshgrid(xval_longlist, yval_longlist, data_longlist, x_meshrgid, y_meshgrid,
                         interp_kind="nearest", fill_val="NaN"):

    return griddata((xval_longlist, yval_longlist),data_longlist,(x_meshrgid, y_meshgrid),
                         method=interp_kind, fill_value=fill_val)

def meshgrid2uniquearrays(xmesh, ymesh, data_meshes):
    """Convert from a meshgrid to unique arrays. The meshes look like:
        xmesh = [[x1, x2, x3, . . . ], [x1, x2, x3, . . . ], . . . ]
        ymesh = [[y1, y1, y1, . . . ], [y2, y2, y2, . . . ], . . . ]
        data_mesh = [[d(x1, y1), d(x2, y1), . . . ], [d(x1, y1), d(x1, y2), . . . ], . . .]"""
    data_unique_arrays = []

    len_x = len(xmesh[0]); len_y = len(xmesh)
    x_array = np.zeros(len_x); y_array = np.zeros(len_y)
    for i in range(0, len_x):
        x_array[i] = xmesh[0,i]

    for i in range(0, len_y):
        y_array[i] = ymesh[i,0]

    # Get the unique data arrays
    for data_mesh in data_meshes:
        data_unique_array = np.zeros((len_x, len_y))
        for x_idx in range(0, len_x):
            for y_idx in range(0, len_y):
                data_unique_array[x_idx, y_idx] = data_mesh[x_idx, y_idx]
        data_unique_arrays.append(data_unique_array)


    return x_array, y_array, data_unique_arrays


def uniquearrays2meshgrids(xarray, yarray, data_arrays, n_xpoints=100, n_ypoints=100,
                       interp_kind="nearest", fill_val="NaN"):
    """ """

    xval_longlist, yval_longlist, data_longlists = unique_arrays2longlists(xarray, yarray, data_arrays)
    x_mesh, y_mesh, interpolated_data = construct_meshgrid(xval_longlist, yval_longlist, data_longlists,
                            n_xpoints=n_xpoints, n_ypoints=n_ypoints, interp_kind=interp_kind, fill_val=fill_val)

    return x_mesh, y_mesh, interpolated_data

def plot_colormap(x, y, data, data_type, n_xpoints=100, n_ypoints=100,
                       interp_kind="nearest", fill_val="NaN", levels=20, colmap="inferno", title="",
                       xlabel="", ylabel=""):

    if data_type == "longlist":
        # Convert to mesh
        xmesh, ymesh, [data_mesh] = construct_meshgrid(x, y, [data],
                                            n_xpoints=n_xpoints, n_ypoints=n_ypoints, interp_kind=interp_kind, fill_val=fill_val)
    elif data_type == "mesh":
        xmesh = x; ymesh = y; data_mesh = data

    elif data_type == "unique":
        xmesh, ymesh, [data_mesh] = uniquearrays2meshgrids(x, y, [data], n_xpoints=n_xpoints,
                                            n_ypoints=n_ypoints, interp_kind=interp_kind, fill_val=fill_val)
    else:
        print("Data type not recognised!")
        return
    # Make plot
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0.1, 0.1, 0.75, 0.75])
    cbax = fig.add_axes([0.9, 0.1, 0.05, 0.75])

    col_contours = ax.contourf(xmesh, ymesh, data_mesh, levels, cmap=colmap) #np.arange(0,1.01,0.01))
    fig.colorbar(col_contours, cax=cbax)

    fig.suptitle(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    plt.show()

    return
