# -*- coding: utf-8 -*-
from aloe import step
from aloe import world

import numpy

from quantarhei import energy_units
from quantarhei import Manager
from quantarhei.core.units import kB_int
from quantarhei import CorrelationFunction
from quantarhei import Molecule
from quantarhei import Aggregate
from quantarhei.qm import RedfieldRateMatrix, RedfieldRelaxationTensor

@step(r'a homodimer with site transition energies (\d+(?:\.\d+)?) "([^"]*)"')
def homodimer_site_energies(self, energy, units):
    print(energy, units)
    world.senergy = float(energy)
    world.h_units = units
    
@step(r'resonance coupling (\d+(?:\.\d+)?) "([^"]*)"')
def resonance_coupling(self, energy, units):
    print(energy, units)
    world.r_coupl = float(energy)
    world.r_units = units        
    
@step(r'I calculate Redfield relaxation rates')
def redfield_relaxation(self):
    print("Redfield rate calculation")

    
    #
    # Correlation function
    #
    params = {"ftype":    world.ctype,
              "reorg":    world.reorg,
              "cortime":  world.ctime,
              "T":        world.temp,
              "matsubara":world.mats}
              
    # FIXME: also time_units, temperature_units
    with energy_units(world.e_units):
        cf = CorrelationFunction(world.ta,params) 

    
    #
    # Homodimer
    #
    with energy_units(world.h_units):
        en = world.senergy
        m1 = Molecule("mol1", [0.0, en])
        m2 = Molecule("mol2", [0.0, en])
#        m3 = Molecule("mol2", [0.0, en])
        
    m1.set_egcf((0,1), cf)
    m2.set_egcf((0,1), cf)
#    m3.set_egcf((0,1), cf) 
    
    agg = Aggregate("Homodimer",maxband=1)
    
    agg.add_Molecule(m1)
    agg.add_Molecule(m2)
#    agg.add_Molecule(m3)
    
#    with energy_units("1/cm"):
#        Hm = m1.get_Hamiltonian()
#        print(Hm)
#        print(m.convert_energy_2_current_u(Hm._data)) 
        
    with energy_units(world.r_units):
        agg.set_resonance_coupling(0,1,world.r_coupl)
#        agg.set_resonance_coupling(1,2,world.r_coupl)
        
    agg.build()
    
    H = agg.get_Hamiltonian()

#    with energy_units("1/cm"):
#        print(H)
#        print(m.convert_energy_2_current_u(H._data))    
    
    sbi = agg.get_SystemBathInteraction()
        
    RRM = RedfieldRateMatrix(H, sbi)
    
    world.K12 = numpy.real(RRM.data[1,2])
    
    
#    #
#    # Homodimer
#    #
#    with energy_units(world.h_units):
#        en = world.senergy
#        m1 = Molecule("mol1", [0.0, en])
#        m2 = Molecule("mol2", [0.0, en])
##        m3 = Molecule("mol2", [0.0, en])
#        
#    m1.set_egcf((0,1), cf)
#    m2.set_egcf((0,1), cf)
##    m3.set_egcf((0,1), cf)
#        
#    agg = Aggregate("Homodimer",maxband=1)
#    
#    agg.add_Molecule(m1)
#    agg.add_Molecule(m2)
##    agg.add_Molecule(m3)
#    
##    with energy_units("1/cm"):
##        Hm = m1.get_Hamiltonian()
##        print(Hm)
##        print(m.convert_energy_2_current_u(Hm._data)) 
#        
#    with energy_units(world.r_units):
#        agg.set_resonance_coupling(0,1,world.r_coupl)
##        agg.set_resonance_coupling(1,2,world.r_coupl)
#        
#    agg.build()
#    
#    H = agg.get_Hamiltonian()
#
##    with energy_units("1/cm"):
##        print(H)
##        print(m.convert_energy_2_current_u(H._data))    
#    
#    sbi = agg.get_SystemBathInteraction()
#    
#    RRT = RedfieldRelaxationTensor(H, sbi)  
#    
#    world.tK12 = numpy.real(RRT.data[1,1,2,2])
    
    
    

@step(r'I get Redfield relaxation rates from file ([^"]*) with rtol (\d+(?:\.\d+)?)')
def compare_rates(self, filename, rtols):
    pass


@step(r'I get Redfield relaxation rates corresponding to analytical results for a homodimer with rtol (\d+(?:\.\d+)?)')
def compare_rates_analytical(self, rtols):
    
    print("Comparison of redfield rates with analytical results")
    m = Manager()
    
    rtol = float(rtols)

    with energy_units(world.r_units):
        J = m.convert_energy_2_internal_u(world.r_coupl)
    with energy_units(world.e_units):
        lamb = m.convert_energy_2_internal_u(world.reorg) 
    tauc = world.ctime
    kBT = kB_int*world.temp
    
    #print(world.temp, lamb, tauc, J)    
    
    K12 = ((lamb*J/2.0)
            *(1.0 + 1.0/numpy.tanh(J/kBT))/((J**2)*tauc + 1.0/(4.0*tauc)))
    
    print(K12, 1.0/K12)
    print(world.K12, 1.0/world.K12)
    
    print(K12/world.K12)
    
#    print(world.tK12/world.K12)
    
    numpy.testing.assert_allclose(K12, world.K12 ,rtol=rtol) #,atol=atol)    