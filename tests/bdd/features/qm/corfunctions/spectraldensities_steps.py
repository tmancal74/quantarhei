# -*- coding: utf-8 -*-

from aloe import step
from aloe import world

import numpy

from quantarhei import energy_units
from quantarhei import Manager
from quantarhei import CorrelationFunction
from ...stepslib import read_n_columns


@step(r'correlation function parameters')
def correlation_function_parameters(self):
    
    for row in self.hashes:
        ftype = row['cf_type']
        temp = float(row['temp'])
        T_units = row['T_units']
        reorg = float(row['reorg'])
        e_units = row['e_units']
        ctime = float(row['ctime'])
        t_units = row['t_units']  
        mats = int(row['mats'])
        
    #params = dict(ftype=ftype,T=temp,reorg=reorg,cortime=ctime,matsubara=mats)
    
    world.e_units = e_units
    world.temp_units = T_units
    world.time_units = t_units
    #world.params = params
    world.ctype = ftype
    world.reorg = reorg
    world.ctime = ctime
    world.temp = temp
    world.mats = mats
    
    
    
@step(r'spectral density is created from correlation function')
def spectral_dens_from_corrfce(self):
    #params = world.params
    
    params = {"ftype":    world.ctype,
              "reorg":    world.reorg,
              "cortime":  world.ctime,
              "T":        world.temp,
              "matsubara":world.mats}
    ta = world.ta
    
    with energy_units(world.e_units):
        cf = CorrelationFunction(ta,params)
        
    world.sd = cf.get_SpectralDensity()

@step(r'spectral density corresponds to file ([^"]*) in internal units')
def compare_data_with_file(self, file):

    print("comparing with file ", file)
    sd_data = read_n_columns(__package__,file,2)
    i = 0
    data = numpy.zeros((world.sd.axis.data.shape[0],2))
    for t in world.sd.axis.data:
        data[i,0] = t
        data[i,1] = numpy.real(world.sd.data[i])
        #data[i,2] = numpy.imag(world.cf.data[i])
        i += 1
    numpy.testing.assert_allclose(sd_data,data,rtol=1.0e-7)
    
@step(r'spectral density corresponds to analytical result for ([^"]*) in internal units')
def compare_spectral_dens_to_analytical(self, fctype):
    m = Manager()
    i = 0
    sd_data = numpy.zeros((world.sd.axis.data.shape[0],2))
    wa = world.ta.get_FrequencyAxis()
    with energy_units("int"):
        sd_data[:,0] = wa.data
        omega = wa.data
    with energy_units(world.e_units):
        lamb = m.convert_energy_2_internal_u(world.reorg)
    ctime = world.ctime
        
    if fctype == "OverdampedBrownian":
        # Analytical for for the overdamped Brownian spectral density
        sd_data[:,1] = (2.0*lamb/ctime)*omega/(omega**2 + (1.0/ctime)**2)
    else:
        raise Exception()
    
    data = numpy.zeros((world.sd.axis.data.shape[0],2))
    for t in world.sd.axis.data:
        data[i,0] = t
        data[i,1] = numpy.real(world.sd.data[i])
        #data[i,2] = numpy.imag(world.cf.data[i])
        i += 1
    numpy.testing.assert_allclose(sd_data,data,rtol=1.0e-7)