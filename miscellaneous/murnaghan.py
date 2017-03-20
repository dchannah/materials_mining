import sys
from pylab import *  # Note - also gives numpy as np
from scipy.optimize import leastsq


def read_e_v(filename):
    """
    Reads in ASCII volume vs. energy data using Numpy's loadtxt function.
    :param filename: ASCII data in some type of xy format.
    :return: [array of volumes], [array of energies]
    """
    with open(filename, 'r') as datfile:
        data = np.loadtxt(datfile)
    v = data[:, 0]  # The first column of data
    e = data[:, 1]  # The second column of data
    return v, e


def parabola_fit(volumes, energies):
    """
    Fits a parabola to the data - used as a starting guess for the EOS.
    :param volumes: List of volumes
    :param energies: List of energies in the same order as the volumes
    :return: Parabolic fit parameters
    """
    a, b, c = np.polyfit(volumes, energies, 2)  # 2nd order is parabola
    return a, b, c


def build_murnaghan_guess(a, b, c):
    """
    Builds an initial guess for Murnaghan EOS from parabolic fitting params.
    :param a: a from y = ax^2 + bx + c fit
    :param b: b from y = ax^2 + bx + c fit
    :param c: c from y = ax^2 + bx + c fit
    :return: Murnaghan params: v0, e0, b0, bP
    """
    v0 = -b/(2 * a)
    e0 = a*v0**2 + b*v0 + c
    b0 = 2*a*v0
    bP = 4  # Empirically, this is always a good guess for most compounds.
    return e0, b0, bP, v0


def murn_fit(params, vols):
    """
    Given a vector of parameters and volumes, return a vector of energies.
    :param params: List of B-M parameters (v0, e0, b0, bP)
    :param vols: Volumes for which to compute energies.
    :return: List of energies for those volumes based on EOS
    """
    e0 = params[0]
    b0 = params[1]
    bP = params[2]
    v0 = params[3]
    energ = e0 + b0*vols/bP*(((v0/vols)**bP)/(bP-1)+1) - b0*b0/(bP-1.)
    return energ


def objective(pars, energies, vols):
    """
    Function to minimize in the least-squares fit.
    :param pars: EOS parameters
    :param energies: List of energies
    :param vols: List of volumes in the same order
    :return: An error function to be minimized in scipy
    """
    err = energies - murn_fit(pars, vols)
    return err


f = sys.argv[1]  # Read in data
v, e = read_e_v(f)
guess_a, guess_b, guess_c = parabola_fit(v, e)  # Initial parabola guess
murn_guess = build_murnaghan_guess(guess_a, guess_b, guess_c)
fitted_params, ier = leastsq(objective, murn_guess, args=(e, v))
print "Bulk modulus " + str(fitted_params[1]*160.21773)
print "First derivative of bulk modulus " + str(fitted_params[2])