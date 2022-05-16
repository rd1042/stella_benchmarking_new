"""There are several different ways to calculate Omega, and
the 'best way' is undetermined."""

from helper_ncdf import extract_data_from_ncdf
import numpy as np
from scipy.optimize import curve_fit
import math
import sys
#from extract_sim_data import get_phiz_data_gs2, get_aparz_data_gs2, get_bparz_data_gs2

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


def calculate_omega_gs2_from_outnc(sim_name):
    """Calculates the growth rate for a simulation bsaed on phi**2(t) and
    also finds the values calculated by GS2."""

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
    except Exception:
        return np.NaN, np.NaN, np.NaN, np.NaN
    fitted_phi2, fitted_growth_rate, growth_rate_error = (
                            fit_growth_rate_from_phi2(t_chopped, phi2_chopped))
    return t, omega_average[:,0,0], fitted_growth_rate, growth_rate_error

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
    if nsteps < MIN_NSTEPS:
        print("Non-critical warning! Number of phi2 data points = ", nsteps)
        print("t = ", t)
        print("phi2 = ", phi2)
    if nsteps < MIN_NSTEPS_ERROR:
        print("Critical warning! Number of phi2 data points = ", nsteps)
        print("t = ", t)
        print("phi2 = ", phi2)
        raise Exception

    return t[(nsteps//2):], phi2[(nsteps//2):] # // because we want an integer, not a float

def write_gs2_data_to_plaintext(sim_longname):
    """For a GS2 sim, write:
    (1) omega(t)
    (2) calcalated omega    # NB this is calculated, but currently not written anywhere
    (3) fields data in stella units """
    outnc_longname = sim_longname   + ".out.nc"

    # Get omega values
    t, omega_average, fitted_growth_rate, growth_rate_error = calculate_omega_gs2_from_outnc(outnc_longname)

    # Write omega(t) and final omega to file
    omega_file_longname = sim_longname + ".omega"
    omega_file = open(omega_file_longname, "w")
    omega_file.write("time \t \t Re[omavg] \t \t Im[omavg] \n")
    if np.isnan(t).any():
        omega_file.write("Error! NaN values detected - probably too few phi2(t) points")
    else:
        for t_idx in range(0, len(t)):
            omega_line = "{:.6e} \t \t {:.12e} \t \t {:.12e} \n".format(t[t_idx], omega_average[t_idx].real, omega_average[t_idx].imag)
            #omega_line = str(t[t_idx]) + "\t \t" + str(omega_average[t_idx].real) + "\t \t" + str(omega_average[t_idx].imag) + "\n"
            omega_file.write(omega_line)

    omega_file.close()
    # Get the fields data
    ## Actually, this is already written to .fields, so don't need to write it.
    theta, phi = extract_data_from_ncdf((sim_longname + ".out.nc"), "theta","phi")
    phi = phi[0][0]
    # print("theta, phi = ", theta, phi)
    # sys.exit()
    try:
        [apar] = extract_data_from_ncdf((sim_longname + ".out.nc"), "apar")
        apar = apar[0][0]/2
    except KeyError:
        apar = np.zeros(len(theta))
    try:
        [bpar] = extract_data_from_ncdf((sim_longname + ".out.nc"), "bpar", "bmag")
        bpar = bpar[0][0]*bmag
    except KeyError:
        bpar = np.zeros(len(theta))

    fields_file_longname = sim_longname + ".final_fields_stella_normalisation"
    fields_file = open(fields_file_longname, "w")
    fields_file.write("theta \t \t Re[phi] \t \t Im[phi] \t \t " +
            "Re[apar] \t \t Im[apar] \t \t Re[apar] \t \t Im[bpar] \t \t Re[bpar] \n")
    for theta_idx in range(0, len(theta)):
        fields_line = "{:.6e} \t \t {:.12e} \t \t {:.12e} \t \t {:.12e} \t \t {:.12e} \t \t {:.12e} \t \t {:.12e} ".format(
                theta[theta_idx], phi[theta_idx].real, phi[theta_idx].imag,
                apar[theta_idx].real, apar[theta_idx].imag,
                bpar[theta_idx].real, bpar[theta_idx].imag)
        fields_file.write(fields_line)
    fields_file.close()

    return
