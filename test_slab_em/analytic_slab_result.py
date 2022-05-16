""" """

from scipy.special import wofz, iv
from scipy.optimize import fsolve, brentq, newton
import numpy as np
import matplotlib.pyplot as plt
from mpmath import findroot

def Gammafunc(arg):
    """Returns e^x * I_0(x), which is the result of integrating J0^2 over vperp"""
    return np.exp(-arg)*iv(0,arg)

def plasmadispfunc(arg, with_diags=True):
    """Returns the plasma dispersion function Z(x) """
    # print("arg, type = ", arg, type(arg))
    # print("wofz(arg):")
    # print(wofz(arg))
    if with_diags:
        print("arg, wofz=", arg, wofz(complex(arg)))
    return wofz(complex(arg))*1j*np.pi#*np.sqrt(np.pi)

def test_plasmadispfunc():
    """ """
    print("plasmadispfunc(1e-5) = ", plasmadispfunc(1e-5))
    print("plasmadispfunc(0.01) = ", plasmadispfunc(0.01))
    print("plasmadispfunc(1) = ", plasmadispfunc(1))
    print("plasmadispfunc(100) = ", plasmadispfunc(100))

class SlabDispersionRelation:

    def __init__(self, Ti_norm=None, Te_norm=None, mi_norm=None, me_norm=None,
                       B_norm=None, beta=None, kpa_norm=None, kperp_norm=None):
        """ """
        self.Ti_norm = Ti_norm
        self.Te_norm = Te_norm
        self.mi_norm = mi_norm
        self.me_norm = me_norm
        self.ni_norm = 1
        self.ne_norm = 1
        self.Zi = 1
        self.Ze = -1
        self.B_norm = B_norm
        self.beta = beta
        self.kpa_norm = kpa_norm
        self.kperp_norm = kperp_norm

        ## Now get the derived values
        self.vthi_norm = (self.Ti_norm/self.mi_norm)**0.5 # Thermal velocity for ions. Normalised to vthr = sqrt(2Tr/mr)
        self.vthe_norm = (self.Te_norm/self.me_norm)**0.5 # Thermal velocity for electrons. Normalised to vthr = sqrt(2Tr/mr)
        self.valfven_norm = (self.B_norm**2 /(self.ni_norm*self.mi_norm*self.beta))**0.5  # Normalised Alfven velocity
        self.Omegai_norm = self.B_norm/self.mi_norm # Gyrofrequency for ions. Normalised: Omegai = Omegai_norm*e*Br/(mr*c)
        self.Omegae_norm = self.B_norm/self.me_norm # Gyrofrequency for electrons. Normalised: Omegai = Omegai_norm*e*Br/(mr*c)
        self.rhoi_norm = self.vthi_norm/self.Omegai_norm # Ion Larmor radius. Normalised: rhoi = rhoi_norm * (vthr/Omegar)
        self.rhoe_norm = self.vthe_norm/self.Omegae_norm # Electron normalised Larmor radius
        self.bi_norm = self.kperp_norm*self.kperp_norm*self.rhoi_norm*self.rhoi_norm/2 # Arugment of the Gamma function for ions
        self.be_norm = self.kperp_norm*self.kperp_norm*self.rhoe_norm*self.rhoe_norm/2 # Arugment of the Gamma function for ions
        self.Gammai = Gammafunc(self.bi_norm)
        self.Gammae = Gammafunc(self.be_norm)
        return

    def dispersion_relation(self, omega_norm, with_diags=True):
        """The LHS of the dispersion relation (RHS=0) in a slab. Derived
        by """
        ### Calculate some normalised quantities
        omegabari = omega_norm/(self.kpa_norm*self.vthi_norm)
        omegabare = omega_norm/(self.kpa_norm*self.vthe_norm)
        dispi = plasmadispfunc(omegabari, with_diags=with_diags)
        dispe = plasmadispfunc(omegabare, with_diags=with_diags)
        one_plus_ombardispi = 1 + omegabari*dispi
        one_plus_ombardispe = 1 + omegabare*dispe

        A0_sum = (self.Zi*self.Zi*self.ni_norm/self.Ti_norm * self.Gammai*one_plus_ombardispi +
                  self.Ze*self.Ze*self.ne_norm/self.Te_norm * self.Gammae*one_plus_ombardispe)
        A0_numeratori = -2*self.beta*omega_norm/(self.kpa_norm*self.vthi_norm)*A0_sum
        A0_numeratore = -2*self.beta*omega_norm/(self.kpa_norm*self.vthe_norm)*A0_sum
        A0_denom = self.kperp_norm*self.kperp_norm - 2*self.beta*(omega_norm/(self.kpa_norm*self.B_norm))**2*A0_sum
        ion_terms = self.Zi*self.Zi*self.ni_norm/self.Ti_norm * (self.Gammai*omegabari*(dispi - (A0_numeratori/A0_denom)*one_plus_ombardispi)+1)
        electron_terms = self.Ze*self.Ze*self.ne_norm/self.Te_norm * (self.Gammae*omegabare*(dispe - (A0_numeratore/A0_denom)*one_plus_ombardispe)+1)
        # print("omega_norm, ion_terms + electron_terms = ", omega_norm, (ion_terms + electron_terms))
        return (ion_terms + electron_terms)

    # def test_func(self, omega_norm):
    #     """ """
    #
    #     omega_alfven = 2
    #     return omega_norm - omega_alfven

    def scan_omega(self, nguesses, low_exponent, high_exponent):
        """A scan of omega, with im(omega)=0 everywhere """
        omega_saw = self.valfven_norm*self.kpa_norm
        guess_array = omega_saw*np.logspace(low_exponent, high_exponent)
        residue_array = np.zeros((nguesses), dtype="complex")
        dispi_array = np.zeros((nguesses), dtype="complex")
        dispe_array = np.zeros((nguesses), dtype="complex")
        one_plus_ombardispi_array = np.zeros((nguesses), dtype="complex")
        one_plus_ombardispe_array = np.zeros((nguesses), dtype="complex")
        for omega_guess_idx in range(0, len(guess_array)):
            omega = guess_array[omega_guess_idx]
            omegabari = omega/(self.kpa_norm*self.vthi_norm)
            omegabare = omega/(self.kpa_norm*self.vthe_norm)
            residue_array[omega_guess_idx] = self.dispersion_relation(guess_array[omega_guess_idx], with_diags=False)
            dispi_array[omega_guess_idx] = plasmadispfunc(omegabari, with_diags=False)
            dispe_array[omega_guess_idx] = plasmadispfunc(omegabare, with_diags=False)
            one_plus_ombardispi_array[omega_guess_idx] = 1+omegabari*plasmadispfunc(omegabari, with_diags=False)
            one_plus_ombardispe_array[omega_guess_idx] = 1+omegabare*plasmadispfunc(omegabare, with_diags=False)

        print("residue_array = ", residue_array)
        print("dispi_array = ", dispi_array)
        print("dispe_array = ", dispe_array)
        print("one_plus_ombardispi_array = ", one_plus_ombardispi_array)
        print("one_plus_ombardispe_array = ", one_plus_ombardispe_array)
        return


if __name__ =="__main__":
    print("Hello world")
    ## Create a "cold ion, long-wavelength" sim.
    my_disp = SlabDispersionRelation(Ti_norm=0.0001, Te_norm=1, mi_norm=1, me_norm=2.7e-4,
                                     B_norm=1, beta=0.1, kpa_norm=1, kperp_norm=0.001)
    ## Expect the sim to converge to shear alfven wave for kperp -> 0
    omega_saw = my_disp.valfven_norm*my_disp.kpa_norm
    result = my_disp.dispersion_relation(omega_saw, with_diags=False)
    print("LHS for omega=SAW omega = ", result)

    ## Root-finding using Newton's method. Fails with NaNs.
    # try:
    #     new_result = newton(my_disp.dispersion_relation, omega_saw, maxiter=50, tol=0.01)
    #     print("new_result = ", new_result)
    #     print("omega, residue = ", new_result, my_disp.dispersion_relation(new_result, with_diags=False))
    # except RuntimeError: # Occurs when newton fails to converge
    #     print("Failed to find root")
    #new_result = findroot(my_disp.dispersion_relation, omega_saw)


    ## Try scanning real values of omega
    # my_disp.scan_omega(100, -2, 2)
