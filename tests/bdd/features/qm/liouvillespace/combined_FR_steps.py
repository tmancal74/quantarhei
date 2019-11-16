# -*- coding: utf-8 -*-

"""

    Test of the combinded standard Foerster and Redfield equations
    in secular approximation.
    


"""
from aloe import step
from aloe import world

import numpy

from quantarhei.testing.feature import FeatureFileGenerator
from quantarhei.testing.feature import match_number 

from quantarhei import Molecule
from quantarhei import Aggregate
from quantarhei import CorrelationFunction
from quantarhei import TimeAxis
from quantarhei import ReducedDensityMatrix
from quantarhei import energy_units
from quantarhei import eigenbasis_of
from quantarhei import Hamiltonian

feature = """
#
#  Combined Redfield and Foerster relaxation theories 
#
#
#

Feature: Combined Redfield-Foerster relaxation theory 

    Combination of Redfield and Foerster relaxation tensors 

@redfield-foerster
#@in_development
Scenario Outline: Redfield-Foerster rates for a small chain of molecules

"""

example_1 = """
    Examples:
       | temp | time_step | nsteps | matsu  | cutoff_coup |atol   |
       | 300  | 0.1       | 10000  | 100    |   50        | 0.01  |
       | 200  | 0.1       | 10000  | 100    |   50        | 0.02  |
       | 100  | 0.1       | 10000  | 100    |   50        | 0.02  |
       | 50   | 0.1       | 10000  | 100    |   50        | 0.02  |

"""


@step(r'time interval from zero with time step '+match_number+r' fs and '
      +match_number+r' steps')
def time_interval(self, time_step, nsteps):
    """Matching the following
    
    #begin_feature    
    Given time interval from zero with time step <time_step> fs and <nsteps> steps
    #end_feature
    
    """
    world.time = TimeAxis(0.0, int(nsteps), float(time_step))
    
def homochain_aggregate_1_redfoer(self):
    """Matching the following 
    
    #begin_feature
    And a small homochain of N molecules with parameters
       | tr_energy | neighbor_coupling | corfce | reorg | ctime | matsubara | temperature | e_units |
       | 12000     |  100              | OB     | 40    | 100   |  <matsu>  |  <temp>     | 1/cm    |
       | 12000     |  80               | OB     | 40    | 100   |  <matsu>  |  <temp>     | 1/cm    | 
       | 12200     |  30               | OB     | 40    | 80    |  <matsu>  |  <temp>     | 1/cm    |
       | 12100     |  200              | OB     | 50    | 100   |  <matsu>  |  <temp>     | 1/cm    |
       | 12100     |  100              | OB     | 50    | 100   |  <matsu>  |  <temp>     | 1/cm    |
    #end_feature       
    """
    pass

@step(r'a small homochain of N molecules with parameters')
def homochain_aggregate(self):    
    molecules = []
    time = world.time #TimeAxis(0.0, 10000, 1.0)
    k = 0
    couplings = []
    temp = 0.0
    for row in self.hashes:
        tenergy = float(row['tr_energy'])
        coupling = float(row['neighbor_coupling'])
        ftype_abbv = row['corfce']
        if ftype_abbv == "OB":
            ftype = "OverdampedBrownian"
        reorg = float(row['reorg'])
        e_units = row['e_units']
        ctime = float(row['ctime'])
        tempp = float(row['temperature'])
        matsu = int(row['matsubara'])
        if temp == 0.0:
            temp = tempp
        else:
            if temp != tempp:
                raise Exception("Temperatures must be the same for all molecules")
        
        params = dict(ftype=ftype, cortime=ctime,
                      reorg=reorg, T=temp, matsubara=matsu)
        
        with energy_units(e_units):
            m = Molecule([0.0, tenergy])
            cf = CorrelationFunction(time, params)
            m.set_transition_environment((0,1), cf)
            
        molecules.append(m)
        couplings.append(coupling)
        k += 1
        
    agg = Aggregate(molecules=molecules)
    
    with energy_units("1/cm"):
        for i in range(k):
            if i+1 < k:
                agg.set_resonance_coupling(i,i+1,couplings[i])
            else:
                agg.set_resonance_coupling(i,0,couplings[i])
            
    agg.build()
    
    world.temperature = temp
    world.aggregate = agg
    world.time = time
    world.N = k+1
    
    
@step(r'I calculate aggregate Redfield-Foerster relaxation tensor with cutoff coupling of '+match_number+r' 1/cm')
def combined_redfield_foerster_tensor_agg(self, ccutoff):
    """Creates aggregate Redfield tensor
    
    #begin_feature
    When I calculate aggregate Redfield-Foerster relaxation tensor with cutoff coupling of <cutoff_coup> 1/cm
    #end_feature
    
    """
    agg = world.aggregate
    time = world.time
    
    world.prop = \
    agg.get_ReducedDensityMatrixPropagator(time,
                                           relaxation_theory="standard_Foerster",
                                           secular_relaxation=True)
    world.ccutoff = float(ccutoff)


@step(r'long time Redfield-Foerster populations will correspond to canonical equilibrium with temperature '
      +match_number+r' K with atol '+match_number)
def thermal_population_comparison(self, temp, atol):
    """Creates canonical populations and checks them against the results from relaxation tensor
    
    #begin_feature
    Then long time Redfield-Foerster populations will correspond to canonical equilibrium with temperature <temp> K with atol <atol>
    #end_feature
    
    """
    
    print("Testing Redfield-Foerster theory")
    
    ccutoff = world.ccutoff
    agg = world.aggregate
    rho0 = ReducedDensityMatrix(dim=world.N)
    
    H = world.aggregate.get_Hamiltonian()
    H1 = Hamiltonian(data=H.data.copy())
    
    with energy_units("1/cm"):
        print(H1)
        print(ccutoff, "1/cm", H1.dim)
        H1.remove_cutoff_coupling(ccutoff)
        print(H1)
    
    print("\nExpected thermal population at ", temp, "K")
    with eigenbasis_of(H1):
        rhoT = agg.get_DensityMatrix(condition_type="thermal_excited_state",
                                 relaxation_theory_limit="strong_coupling",
                                 temperature=world.temperature)
        pop_T = numpy.zeros(world.N, dtype=numpy.float64)
    
        for i in range(world.N):
            pop_T[i] = numpy.real(rhoT.data[i,i])
            print(i, ":", pop_T[i])
            
    print("Final calculated population")
    rho0.data[world.N-1,world.N-1] = 1.0
    rho_t = world.prop.propagate(rho0)
    pop_C = numpy.zeros(world.N, dtype=numpy.float64)
    for i in range(world.N):
        pop_C[i] = numpy.real(rho_t.data[world.time.length-1,i,i])
        print(i, ":", pop_C[i])
    
    print("absolute tolerance atol=", float(atol))
    numpy.testing.assert_allclose(pop_C, pop_T, atol=float(atol))


                                           
if __name__ == '__main__':
    """
    
    Feature file: combined_FR_1.feature
    
    """
    gen = FeatureFileGenerator(feature+example_1)
    
    gen.add_Given(time_interval)
    gen.add_Given(homochain_aggregate_1_redfoer)
    gen.add_When(combined_redfield_foerster_tensor_agg)
    gen.add_Then(thermal_population_comparison)
    
    gen.generate_feature_file("combined_FR_1.feature")
                                         