# -*- coding: utf-8 -*-
"""

    Demonstration of Quantarhei fitting capabilities


"""

_show_plots_ = True

import numpy
import quantarhei as qr
import matplotlib.pyplot as plt

###############################################################################
#
#           Fitting by exponential functions
#
###############################################################################

#
# Example functions to fit
#

def single_exp(x):
    """Example single exponential function
    
    """
    
    return numpy.exp(-x/245.0) - 0.2


def double_exp(x):
    """Example double exponential function
    
    
    """
    
    return 0.5*numpy.exp(-x/35.0) + 0.4*numpy.exp(-x/801.0) + 0.1


#
#  Reference an test functions
#

def ref_nexp(x, par):
    """This is how exponential function is constructed from fitting parameters
    
    
    """
    y = numpy.zeros(x.shape[0], dtype=x.dtype)
    for ii in range(0, len(par)-1, 2):
        y += par[ii]*numpy.exp(-par[ii+1]*x)
    y += par[ii+2]
    return y

    
def test_exponential(x, f, inpar, fitpar):
    """
    
    """
    y = ref_nexp(x, inpar)
    if _show_plots_:
        plt.plot(x, y, "-r")
        f.plot(show=False)
        plt.plot(x, ref_nexp(x, fitpar), "--g")
        plt.show()


#
# Time axis on which we will fit
#
t = qr.TimeAxis(0.0, 1000, 1.0)
x = t.data


#
# Mono-exponential function
#
fv = single_exp(x)
f = qr.DFunction(t, fv)

# parameters a, b and c of a*exp(-b*t) + c

inpar = [1.5, 1.0/100.0, 0.1]

# This is how we use the `fit_exponential` function 
fitpar = f.fit_exponential(guess=inpar)

test_exponential(x, f, inpar, fitpar)

#
# Bi-exponential function
#
fv = double_exp(x)
f = qr.DFunction(t, fv)

# pairs of params a_n, b_n and one parameter c of \sum_n a_n*exp(-b_n*t) + c
inpar = [1.5, 1.0/100.0, 1.0, 1.0/500, 0.1]

# This is how we use the `fit_exponential` function
fitpar = f.fit_exponential(guess=inpar)

test_exponential(x, f, inpar, fitpar)


###############################################################################
#
#           Fitting by Gaussain functions
#
###############################################################################

#
# Example functions to fit
#

def single_gauss(x):
    """Example single exponential function
    
    """
    
    return 1.34*numpy.exp(-(x-55.3)**2/(24.0**2)) - 0.2


def double_gauss(x):
    """Example double exponential function
    
    
    """
    
    return (0.5*numpy.exp(-(x-34.6)**2/(10.0**2)) + 
            0.4*numpy.exp(-(x-71.0)**2/(18.0**2)) + 0.1)


#
#  Reference an test functions
#

def ref_ngauss(x, par):
    """This is how exponential function is constructed from fitting parameters
    
    
    """
    y = numpy.zeros(x.shape[0], dtype=x.dtype)
    for ii in range(0, len(par)-1, 3):
        y += par[ii]*numpy.exp(-4.0*numpy.log(2.0)*(x-par[ii+1])**2
                               /(par[ii+2]**2))
    y += par[ii+3]
    return y


def test_gaussian(x, f, inpar, fitpar):
    """
    
    """
    y = ref_ngauss(x, inpar)
    if _show_plots_:
        plt.plot(x, y, "-r")
        f.plot(color="-k", show=False)
        plt.plot(x, ref_ngauss(x, fitpar), "--g")
        plt.show()

#
# Time axis on which we will fit
#
t = qr.TimeAxis(0.0, 100, 1.0)
x = t.data

#
# Single Gaussian function
#
fv = single_gauss(x)
f = qr.DFunction(t, fv)
inpar = [1.5, 50.0, 100.0, 0.1]

# This is how we use the `fit_exponential` function 
fitpar = f.fit_gaussian(guess=inpar)

test_gaussian(x, f, inpar, fitpar)

#
# Two Gaussian functions
#
fv = double_gauss(x)
#f = qr.DFunction(t, fv)
f = qr.PumpProbeSpectrum()
f.set_axis(t)
f.set_data(fv)

inpar = [1.5, 30.0, 20.0, 0.6, 60.0, 20.0, 0.1]

# This is how we use the `fit_exponential` function 
fitpar = f.fit_gaussian(guess=inpar)

test_gaussian(x, f, inpar, fitpar)