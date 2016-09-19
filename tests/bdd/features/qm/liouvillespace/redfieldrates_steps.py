# -*- coding: utf-8 -*-
from aloe import step
from aloe import world

from ...stepslib import upper_half_time_axis

@step(r'a homodimer with site transition energies (\d+(?:\.\d+)?) "([^"]*)"')
def homodimer_site_energies(self, energy, units):
    print(energy, units)
    world.senergy = float(energy)
    world.e_units = units
    
@step(r'resonance coupling (\d+(?:\.\d+)?) "([^"]*)"')
def resonance_coupling(self, energy, units):
    print(energy, units)
    world.r_coupl = float(energy)
    world.r_units = units
    
@step(r'correlation function with parameters')
def correlation_function(self):
    for row in self.hashes:
        ctype = row['cf_type']
        
    
@step(r'I calculate Redfield relaxation rates')
def redfield_relaxation(self):
    pass

@step(r'I get rates from file ([^"]*)')
def compare_rates(self,filename):
    pass
