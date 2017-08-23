# -*- coding: utf-8 -*-

from aloe import step
from aloe import world
import numpy
import pkg_resources
from unittest import TestCase

from quantarhei import TimeAxis
from quantarhei import FrequencyAxis


@step(r'upper-half TimeAxis with parameters')
def upper_half_time_axis(self):
    
    for row in self.hashes:
        start = float(row['start'])
        Ns    = int(row['number_of_steps'])
        dt    = float(row['step'])
        
    print("TimeAxis", start, Ns, dt)
    # upper-half is default
    world.ta = TimeAxis(start,Ns,dt)
    
    tc = TestCase()
    tc.assertEqual(world.ta.atype,"upper-half")

@step(r'complete TimeAxis with parameters')
def complete_time_axis(self):
    
    for row in self.hashes:
        start = float(row['start'])
        Ns    = int(row['number_of_steps'])
        dt    = float(row['step'])
        
    print("TimeAxis", start, Ns, dt, "atype='complete'")
    # upper-half is default
    world.ta = TimeAxis(start,Ns,dt,atype="complete")
    
    tc = TestCase()
    tc.assertEqual(world.ta.atype,"complete")    
 
@step(r'complete FrequencyAxis with parameters')
def complete_frequency_axis(self):
    
    for row in self.hashes:
        start = float(row['start'])
        Ns    = int(row['number_of_steps'])
        dw    = float(row['step'])
        
    print("FrequencyAxis", start, Ns, dw,"atype='complete'")
    world.wa = FrequencyAxis(start,Ns,dw)
    
    tc = TestCase()
    tc.assertEqual(world.wa.atype,"complete")

@step(r'upper-half FrequencyAxis with parameters')
def upper_half_frequency_axis(self):
    
    for row in self.hashes:
        start = float(row['start'])
        Ns    = int(row['number_of_steps'])
        dw    = float(row['step'])
        
    print("FrequencyAxis", start, Ns, dw, "atype=upper-half")
    world.wa = FrequencyAxis(start,Ns,dw,atype="upper-half")

    tc = TestCase()
    tc.assertEqual(world.wa.atype,"upper-half")
    
@step(r'CorrelationFunction with parameters')
def correlation_function(self):
    for row in self.hashes:
        pass
    
    
@step(r'temperature (\d+(?:\.\d+)?) "([^"]*)"')
def temperature(self, temp, T_units):
    print("temperature ", temp, T_units)
    world.temp = float(temp)
    world.T_units = T_units    
    
def read_n_columns(path, file, n):
    """Reads n columns of data from a package file """
    data = pkg_resources.resource_string(path,file)    
    A = numpy.fromstring(data.decode('utf-8'),dtype=numpy.float32,sep=' ')
    N = A.shape[0]
    Nh = int(N/n)
    try:
        ret = A.reshape(Nh,n)
    except:
        raise Exception()
        
    return ret    
