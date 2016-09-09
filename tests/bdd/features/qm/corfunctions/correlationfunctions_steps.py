# -*- coding: utf-8 -*-
from aloe import step
from aloe import world

import pkg_resources
import numpy

from quantarhei import TimeAxis
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

@step(r'temperature (\d+(?:\.\d+)?) "([^"]*)"')
def temperature(self, temp, T_units):
    print("temperature ", temp, T_units)
    world.temp = float(temp)
    world.T_units = T_units

@step(r'number of Matsubara frequencies (\d+(?:\.\d+)?)')
def matsubara(self, Nm):
    print("no. Matsubara frequencies ", Nm)
    world.mats = int(Nm)


@step(r'time interval')
def time_axis(self):
    
    for row in self.hashes:
        start = float(row['start'])
        Ns    = int(row['number_of_steps'])
        dt    = float(row['step'])
        
    print("TimeAxis", start, Ns, dt)
    world.ta = TimeAxis(start,Ns,dt)

@step(r'When I calculate the ([^"]*) correlation function')
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
        world.cf = CorrelationFunction(world.ta,params)    
    
    
def read_n_columns(file,n):
    """Reads n columns of data from a package file """
    data = pkg_resources.resource_string(__package__,file)    
    A = numpy.fromstring(data.decode('utf-8'),dtype=numpy.float32,sep=' ')
    N = A.shape[0]
    Nh = int(N/n)
    try:
        ret = A.reshape(Nh,n)
    except:
        raise Exception()
        
    return ret    


@step(r'Then I get data from the file ([^"]*) in internal units')
def compare_data_with_file(self, file):

    print("comparing with file ", file)
    cf_data = read_n_columns(file,2)
    print(cf_data)