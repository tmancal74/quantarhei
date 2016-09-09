# -*- coding: utf-8 -*-

from aloe import step
from aloe import world

from quantarhei import TimeAxis

@step(r'TimeAxis')
def time_axis(self):
    
    for row in self.hashes:
        start = float(row['start'])
        Ns    = int(row['number_of_steps'])
        dt    = float(row['step'])
        
    print("TimeAxis", start, Ns, dt)
    world.ta = TimeAxis(start,Ns,dt)
    
@step(r'CorrelationFunction')
def correlation_function(self):
    for row in self.hashes:
        pass
    
    
    