# -*- coding: utf-8 -*-
from aloe import step
from aloe import world

from ...stepslib import upper_half_time_axis
from ...stepslib import temperature
from ...stepslib import read_n_columns

import numpy

#from quantarhei import TimeAxis
from quantarhei import CorrelationFunction
from quantarhei import energy_units

@step(r'reorganization energy (\d+(?:\.\d+)?) "([^"]*)"') 
def reorganization_energy(self, reorg, e_units):
    print("\nreorganization energy ", reorg, e_units)
    world.reorg = float(reorg)
    world.e_units = e_units

@step(r'correlation time (\d+(?:\.\d+)?) "([^"]*)"')
def correlation_time(self, ctime, t_units):
    print("correlation time ", ctime, t_units)
    world.ctime = float(ctime)
    world.t_units = t_units


@step(r'number of Matsubara frequencies (\d+(?:\.\d+)?)')
def matsubara(self, Nm):
    print("no. Matsubara frequencies ", Nm)
    world.mats = int(Nm)

@step(r'correlation function is created from spectral density')
def corfce_from_spectdens(self):

    sd = world.sd        
    world.cf = sd.get_CorrelationFunction()


#@step(r'spectral density is back-created from correlation function')
#def spectral_dens_back_from_corrfce(self):
#
#    cf = world.cf        
#    world.sd = cf.get_SpectralDensity()

@step(r'I calculate the ([^"]*) correlation function')
def correlation_function_of_type(self, ctype):
    print("correlation function type ", ctype)
    world.ctype = ctype

    params = {"ftype":    world.ctype,
              "reorg":    world.reorg,
              "cortime":  world.ctime,
              "T":        world.temp,
              "matsubara":world.mats}
              
    # FIXME: also time_units, temperature_units
    with energy_units(world.e_units):
        cf = CorrelationFunction(world.ta,params) 
    
    world.cf = cf
    

@step(r'correlation function corresponds to file ([^"]*) in internal units with rtol (\d+(?:\.\d+)?) and atol (\d+(?:\.\d+)?)')
def compare_data_with_file(self, file, rtols, atols):

    rtol = float(rtols)
    atol = float(atols)
    print("comparing with file ", file)
    cf_data = read_n_columns(__package__,file,3)
    i = 0
    data = numpy.zeros((world.ta.data.shape[0],3))
    for t in world.ta.data:
        data[i,0] = t
        data[i,1] = numpy.real(world.cf.data[i])
        data[i,2] = numpy.imag(world.cf.data[i])
        i += 1
    numpy.testing.assert_allclose(cf_data,data,rtol=rtol,atol=atol)