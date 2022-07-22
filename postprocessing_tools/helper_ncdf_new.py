"""Helper module, containing often-used and mature methods. The aim is for this
script to be thoroughly tested, and then called by other scripts as needed.

Author: Bob Davies"""

import ncdf2dict
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import sys
import xarray as xr

LINEAR_SIMS = "../../gs2_sims_linear/"


def linear_growth(xdata, a, b):
    return a + b*xdata

def fit_growth_rate_from_phi2(t, phi2):
    """Fit a straight line to the natural logarithm of phi**2. If we assume that
    phi**2 is described by phi**2 = A*exp(Bt), then log(phi**2) = ln(A) + Bt"""

    logphi2 = np.log(phi2)

    popt, pcov = curve_fit(linear_growth, t, logphi2)
    [a_opt, b_opt] = popt
    growth_rate = b_opt/2
    growth_rate_stdev = (np.sqrt(pcov[1, 1]))/2
    fitted_phi2 = np.exp(a_opt)*np.exp(b_opt*t)

    return fitted_phi2, growth_rate, growth_rate_stdev

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

def view_ncdf_variables(outnc_longname):
    """View the names all variables in the netcdf file"""
    data = ncdf2dict.ncdf2dict(outnc_longname)
    print(list(data.keys()))
    return

def compare_growth_rate_to_gs2(fitted_growth_rate, omega_average):
    """Compares the fitted value of the growth rate to the values given
    by GS2."""

    OMEGA_CONV_TOLERANCE = 10**-4   # Maximum st. dev. of growth rates sample
    OMEGA_FIT_TOLERANCE = 10**-3    # Maximum fractional difference between fitted
                                    # and GS2-provided growth rates.
    warnings = 0
    growth_rates = omega_average.imag
     # We find the converged frequency from the final 20% of the frequency
     # values
    sample_size = max(len(growth_rates)//5, 20) # // to force sample_size to integer value
    if len(growth_rates) < sample_size:
       # print("Warning! Not enough growth rate values to compare with GS2 reliably.")
       sample_size = len(growth_rates) - 1

    sample_growth_rates = growth_rates[-sample_size:]

    converged_growth_rate = np.average(sample_growth_rates)
    sample_stdev = np.std(sample_growth_rates)
    # if sample_stdev > OMEGA_CONV_TOLERANCE:
    #     print(("Warning! Growth may not not converged, standard deviation " +
    #            "frequency sample is " + str(sample_stdev)))

    # Compare converged and fitted growth rates - we find the fractional
    # difference by dividing by the converged growth rate rather than the
    # fitted growth rate, because fitting errors may cause the latter to have
    # strange values.
    # if (((converged_growth_rate - fitted_growth_rate)/converged_growth_rate)
    #     > OMEGA_FIT_TOLERANCE):
    #     warnings +=1
    #     print(("Warning! Fitted and converged growth rates do not match " +
    #            "within specified tolerance. Converged growth rate, " +
    #            "fitted growth rate = " + str(converged_growth_rate) + ", " +
    #            str(fitted_growth_rate)))

    # if warnings == 0:
    #     print("Growth rate fit successful!")
    return

def calculate_converged_frequency(t, omega_average):
    """Calculates the frequency given by GS2, and checks if it has converged to
    a reasonable level.

    Bob 23/04/21: Modified to simply take the last point for frequency; frequency tends to
    bifurcate, rather than converge, so the only way to see if we're actually
    converge is to look at frequency(time)"""

    FREQ_CONV_TOLERANCE = 10**-4

    frequency = omega_average.real
     # We find the converged frequency from the final 20% of the frequency
     # values
    sample_size = max(len(frequency)//5, 20) # // to force sample_size to integer value
    if len(frequency) < sample_size:
       #print("Warning! Not enough growth rate values to compare with GS2 reliably.")
       sample_size = len(frequency) -1
    sample_frequencies = frequency[-sample_size:]

    converged_frequency = np.average(sample_frequencies)
    sample_stdev = np.std(sample_frequencies)
    # if sample_stdev > FREQ_CONV_TOLERANCE:
    #     print(("Warning! Frequency may not not converged, standard deviation " +
    #            "frequency sample is " + str(sample_stdev)))

    conv_freq_tlims = [t[-sample_size], t[-1]]
    return frequency[-1][0][0], sample_stdev, conv_freq_tlims

def chop_fitting_data(t, phi2):
    """Chop t and phi2 for fitting the growth rate, according to the following rules:
     - Remove 'nan's and 'inf's
     - Of the remaining data, take the latter half of it """
    MIN_NSTEPS = 20
    MIN_NSTEPS_ERROR = 4
    # Convert from lists to arrays
    t = np.asarray(t)
    phi2 = np.asarray(phi2)

    # Check that all time values are finite - if it isn't, something strange has
    # occured and we should terminate.
    if not np.isfinite(t).all():
        print("Error! Non-finite t values detected.")
        print("t= ", t)
        sys.exit("Script stopped")


    if not np.isfinite(phi2).all():
        #print("Warning! Non-finite phi2 values detected")
        valid_value_array = np.isfinite(phi2)
        invalid_value_idxs = np.where(valid_value_array == False)
        #print "False idxs = ", invalid_value_idxs
        #print "invalid_value_idxs[0] = ", invalid_value_idxs[0]
        t = t[:invalid_value_idxs[0][0]]; phi2 = phi2[:invalid_value_idxs[0][0]]

    nsteps = len(t)
    # if nsteps < MIN_NSTEPS:
        # print("Non-critical warning! Number of phi2 data points = ", nsteps)
        # print("t = ", t)
        # print("phi2 = ", phi2)
    if nsteps < MIN_NSTEPS_ERROR:
        print("Critical warning! Number of phi2 data points = ", nsteps)
        print("t = ", t)
        print("phi2 = ", phi2)
        raise Exception

    return t[(nsteps//2):], phi2[(nsteps//2):] # // because we want an integer, not a float

def lazy_calculate_omega(sim_name, return_fitted_vals=False):
    """Calculate the growth rate based on phi**2(t), designed for simulations
    where omega_average has not been saved. """

    [t, phi2] = extract_data_from_ncdf(sim_name, 't', 'phi2')

    # Only use the last half of the data to fit the growth rate - hopefully,
    # growh rate has converged at this point.
    print("sim_name = ", sim_name)
    t_chopped, phi2_chopped = chop_fitting_data(t, phi2)
    fitted_phi2, fitted_growth_rate, growth_rate_error = (
                            fit_growth_rate_from_phi2(t_chopped, phi2_chopped))

    if return_fitted_vals==True:
       sys.exit("Method  currently has no implementation of return_fitted_vals, aborting.")
    else:
        return (fitted_growth_rate, growth_rate_error)

def calculate_omega(sim_name, return_fitted_vals=False):
    """Calculates the growth rate for a simulation bsaed on phi**2(t) and
    compares it to the values calculated by GS2."""

    try:
        [t, phi2, omega_average] = extract_data_from_ncdf(sim_name, 't', 'phi2',
                                                        'omega_average')
    except KeyError:
        [t, phi2, omega_average] = extract_data_from_ncdf(sim_name, 't', 'phi2',
                                                        'omegaavg')

    # Only use the last half of the data to fit the growth rate - hopefully,
    # growh rate has converged at this point.

    try:
        t_chopped, phi2_chopped = chop_fitting_data(t, phi2)
        fitted_phi2, fitted_growth_rate, growth_rate_error = (
                                fit_growth_rate_from_phi2(t_chopped, phi2_chopped))
        compare_growth_rate_to_gs2(fitted_growth_rate, omega_average)
        converged_frequency, freq_error, conv_freq_tlims = calculate_converged_frequency(t, omega_average)
    except Exception:   # Likely caused by too few points to plot.
        t_chopped = None
        fitted_phi2 = None
        fitted_growth_rate =  omega_average[-1].imag
        growth_rate_error = None
        converged_frequency = omega_average[-1].real
        freq_error = None
        conv_freq_tlims = None

    if return_fitted_vals==True:
       return (t, t_chopped, phi2, fitted_phi2, fitted_growth_rate, growth_rate_error,
                converged_frequency, freq_error, conv_freq_tlims, omega_average)
    else:
        return (fitted_growth_rate, growth_rate_error, converged_frequency,
                freq_error)

def calculate_omega_kykx(sim_name, return_fitted_vals=False):
    """Calculates the growth rate for a simulation bsaed on phi**2(t) and
    compares it to the values calculated by GS2."""


    #view_ncdf_variables(sim_name)
    try:
        [t, phi2_kykx, omega_average_kykx] = extract_data_from_ncdf(sim_name, 't', 'phi2_by_mode',
                                                        'omega_average')
    except KeyError:
        [t, phi2_kykx, omega_average_kykx] = extract_data_from_ncdf(sim_name, 't', 'phi2_by_mode',
                                                        'omegaavg')

    # print("phi2.shape = ", phi2.shape)
    # print("omega_average.shape = ", omega_average.shape)
    n_ky = phi2_kykx.shape[1]; n_kx = phi2_kykx.shape[2]

    fitted_growth_rate_kykx = np.zeros((n_ky, n_kx))
    growth_rate_error_kykx = np.zeros((n_ky, n_kx))
    converged_frequency_kykx = np.zeros((n_ky, n_kx))
    freq_error_kykx = np.zeros((n_ky, n_kx))

    for ky_idx in range(0, n_ky):
        for kx_idx in range(0, n_kx):

            phi2 = phi2_kykx[:, ky_idx, kx_idx]
            omega_average = omega_average_kykx[:, ky_idx, kx_idx]

            # Only use the last half of the data to fit the growth rate - hopefully,
            # growh rate has converged at this point.

            t_chopped, phi2_chopped = chop_fitting_data(t, phi2)
            fitted_phi2, fitted_growth_rate, growth_rate_error = (
                                    fit_growth_rate_from_phi2(t_chopped, phi2_chopped))
            #print("calculate_omega fitted_growth_rate = ", fitted_growth_rate)
            compare_growth_rate_to_gs2(fitted_growth_rate, omega_average)
            converged_frequency, freq_error, conv_freq_tlims = calculate_converged_frequency(t, omega_average)

            fitted_growth_rate_kykx[ky_idx, kx_idx] = fitted_growth_rate
            growth_rate_error_kykx[ky_idx, kx_idx] = growth_rate_error
            converged_frequency_kykx[ky_idx, kx_idx] = converged_frequency
            freq_error_kykx[ky_idx, kx_idx] = freq_error


    #
    # if return_fitted_vals==True:
    #    return (t, t_chopped, phi2, fitted_phi2, fitted_growth_rate, growth_rate_error,
    #             converged_frequency, freq_error, conv_freq_tlims, omega_average)

    return (fitted_growth_rate_kykx, growth_rate_error_kykx, converged_frequency_kykx,
            freq_error_kykx)

def plot_growth_rates(sim_name, title, save_loc="./", include_phi=False):
    """For visual inspection of omega-related outputs of simulation. Plots
    phi**2(t) and omega_average(t), along with the fitted growth rate."""

    (t, t_chopped, phi2, fitted_phi2, fitted_growth_rate, growth_rate_error, converged_frequency,
        freq_error, conv_freq_tlims, omega_average) = calculate_omega(sim_name, return_fitted_vals=True)
    try:
        [theta, phi_by_mode] = extract_data_from_ncdf(sim_name, 'theta', 'phi')
    except KeyError:
        print("Error getting mode structure")
        return

    absphi = phi_by_mode[0][0]
    absphi = abs(absphi/max(abs(absphi)))

    # Plot phi**2(t)
    fig = plt.figure(1, figsize=(17.0, 9.6))
    fig.suptitle(title, fontsize=24, fontweight='bold')
    ax1 = fig.add_subplot(221)
    ax1.scatter(t, phi2, c="black", marker="+", label="phi**2")
    if fitted_phi2 is not None:
        ax1.plot(t_chopped, fitted_phi2, c='r',
            label="Fitted to log(phi**2)")
    ax1.set_yscale('log')
    ax1.legend(loc="best")
    ax1.set_xlabel("time (normalised GS2 units)")
    ax1.set_ylabel("phi**2 (normalised GS2 units)")

    # Plot frequency(t)
    ax2 = fig.add_subplot(222)
    ax2.scatter(t, omega_average.real, c='black', marker ="x", label="omega_average")
    #print("omega_average.shape, omega_average[-1] = ", omega_average.shape, omega_average[-1])
    ax2.text(np.min(t), np.max(omega_average.real)/2, "freq[-1] = {:.5f}".format(omega_average[-1,0,0].real))
    if conv_freq_tlims is not None:
        ax2.plot(conv_freq_tlims, [converged_frequency, converged_frequency], c="r",
                    label="converged frequency")
    #ax2.set_xlabel("time (normalised GS2 units)")
    ax2.set_ylabel("frequency (normalised GS2 units)")
    ax2.legend(loc='best')
    ax2.tick_params('y', which='both', labelsize=8, direction="in")
    ax2.tick_params('x', which='both', labelsize=8, direction="in", bottom=False)
    ax2.grid(True, which="both", axis="both", linestyle="--")
    ax2.set_axisbelow(True)
    plt.setp(ax2.get_xticklabels(), visible=False)

    # Plot GS2's growth rate(t)
    ax3 = fig.add_subplot(224, sharex=ax2)
    ax3.scatter(t, omega_average.imag, c='black', marker ="x", label="omega_average")
    ax3.text(np.min(t), np.max(omega_average.imag)/2, "freq[-1] = {:.5f}".format(omega_average[-1,0,0].imag))
    if t_chopped is not None:
        ax3.plot([t_chopped[0], t_chopped[-1]], [fitted_growth_rate, fitted_growth_rate],
                 c='r', label="Fitted to log(phi**2)")
    ax3.set_xlabel("time (normalised GS2 units)")
    ax3.legend(loc='best')
    ax3.set_ylabel("growth rate (normalised GS2 units)")
    ax3.set_axisbelow(True)
    ax3.tick_params('both', which='both', labelsize=8, direction="in")
    ax3.grid(True, which="both", axis="both", linestyle="--")
    #ax3.title(title)
    ax4 = fig.add_subplot(223)
    ax4.plot(theta/np.pi, absphi, c="black")
    ax4.set_xlabel(r"$\theta/\pi$")
    ax4.set_ylabel(r"$\vert \phi \vert$")
    ax4.tick_params('both', which='both', labelsize=8, direction="in")
    ax4.grid(True, which="both", axis="both", linestyle="--")
    plt.savefig(save_loc + title + ".png")
    plt.close()
    return

def plot_phi2_from_nc(sim_name, title, save_loc="./"):
    """For visual inspection of omega-related outputs of simulation. Plots
    phi**2(t) and omega_average(t), along with the fitted growth rate."""

    [t, phi2, omega_average] = extract_data_from_ncdf(sim_name, 't', 'phi2',
                                                      'omega_average')

    # Plot phi**2(t)
    fig = plt.figure(1, figsize=(17.0, 9.6))
    fig.suptitle(title, fontsize=24, fontweight='bold')
    ax1 = fig.add_subplot(111)
    ax1.scatter(t, phi2, c="black", marker="+", label="phi**2")
    ax1.set_yscale('log')
    ax1.legend(loc="best")
    ax1.set_xlabel("time (normalised GS2 units)")
    ax1.set_ylabel("phi**2 (normalised GS2 units)")

    plt.savefig(save_loc + title + ".png")
    plt.close()

    return

def plot_multiple_growth_rates(sims, title, save_loc="./", labels=[]):
    """For visual inspection of omega-related outputs of simulation. Plots
    phi**2(t) and omega_average(t), along with the fitted growth rate."""


    if len(labels) == 0:
        labels = sims

    fig = plt.figure(1, figsize=(17.0, 9.6))
    fig.suptitle(title, fontsize=24, fontweight='bold')
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(224, sharex=ax2)

    for i, sim_name in enumerate(sims):
        (t, t_chopped, phi2, fitted_phi2, fitted_growth_rate, growth_rate_error, converged_frequency,
            freq_error, conv_freq_tlims, omega_average) = calculate_omega(sim_name, return_fitted_vals=True)

        # Plot phi**2(t)
        ax1.scatter(t, phi2, marker="+")
        ax1.plot(t_chopped, fitted_phi2, label=labels[i])

        # Plot frequency(t)
        ax2.scatter(t, omega_average.real, marker ="x", label=labels[i])
        ax2.plot(conv_freq_tlims, [converged_frequency, converged_frequency])

        # Plot GS2's growth rate(t)
        ax3.scatter(t, omega_average.imag, marker ="x", label=labels[i])
        ax3.plot([t_chopped[0], t_chopped[-1]], [fitted_growth_rate, fitted_growth_rate])

    ax1.set_yscale('log')
    ax1.legend(loc="best")
    ax1.set_xlabel("time (normalised GS2 units)")
    ax1.set_ylabel("phi**2 (normalised GS2 units)")
    ax2.set_ylabel("frequency (normalised GS2 units)")
    ax2.legend(loc='best')
    ax2.tick_params('y', which='both', labelsize=8, direction="in")
    ax2.tick_params('x', which='both', labelsize=8, direction="in", bottom=False)
    ax2.grid(True, which="both", axis="both", linestyle="--")
    ax2.set_axisbelow(True)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax3.set_xlabel("time (normalised GS2 units)")
    ax3.legend(loc='best')
    ax3.set_ylabel("growth rate (normalised GS2 units)")
    ax3.set_axisbelow(True)
    ax3.tick_params('both', which='both', labelsize=8, direction="in")
    ax3.grid(True, which="both", axis="both", linestyle="--")
    plt.show()
    #plt.savefig(save_loc + title + ".png")
    #plt.close()
    return

def view_ncdf_variables_with_xarray(sim_name):
    """View the names all variables in the netcdf file"""
    data = xr.open_dataset(sim_name)
    print(list(data.keys()))
    return
    
if __name__ == "__main__":
    outnc_longname = LINEAR_SIMS + "lowq0_psinkx_scan_millerparams_1/run_psin_0.15_ky_0.15_kx_0.03/input.out.nc"
    view_ncdf_variables(outnc_longname)
    [kx] = extract_data_from_ncdf(outnc_longname, "kx")
    print("kx = ", kx)
