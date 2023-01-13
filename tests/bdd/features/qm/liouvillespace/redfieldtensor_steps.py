# -*- coding: utf-8 -*-
from aloe import step
from aloe import world

import numpy

from quantarhei import energy_units, eigenbasis_of
from quantarhei import CorrelationFunction
from quantarhei import Molecule
from quantarhei import Aggregate
from quantarhei.qm import RedfieldRelaxationTensor

from quantarhei import Manager

@step(r'I calculate Redfield relaxation tensor')
def redfield_tensor(self):
    print("Redfield tensor calculation")

    
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
        m1 = Molecule([0.0, en], "mol1")
        m2 = Molecule([0.0, en], "mol2")
        
    m1.set_egcf((0,1), cf)
    m2.set_egcf((0,1), cf)
        
    agg = Aggregate(name="Homodimer")
    
    agg.add_Molecule(m1)
    agg.add_Molecule(m2)
    
#    m = Manager()
##
#    with energy_units("1/cm"):
#        Hm = m1.get_Hamiltonian()
#        print(Hm)
#        print(m.convert_energy_2_current_u(Hm._data)) 
        
    with energy_units(world.r_units):
        agg.set_resonance_coupling(0,1,world.r_coupl)
        
    agg.build()
    
    H = agg.get_Hamiltonian()
    world.HH = H

##
#    with energy_units("1/cm"):
#        print(H)
#        print(m.convert_energy_2_current_u(H._data))    
    
    sbi = agg.get_SystemBathInteraction()
    
    H.protect_basis()
    with eigenbasis_of(H):    
        RRT = RedfieldRelaxationTensor(H, sbi)
    
        dim = H.dim
        rates_T = numpy.zeros(dim*dim)
        k = 0

        world.K12 = numpy.real(RRT.data[1,1,2,2])
        for i in range(dim):
            for j in range(dim):
                rates_T[k] = numpy.real(RRT.data[i,i,j,j])
                k += 1
                
    world.rates_T = rates_T
    
    