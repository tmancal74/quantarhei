# -*- coding: utf-8 -*-
#******************************************************************************
#
# Simulation of single molecule spectra of RC ...
# using Quantarhei package, version ....
#
# Author: Tomas Mancal
#
# Last changes: July 14th, 2016
#
#
#
#******************************************************************************
#
# About the code
# --------------
#
# The code is explained in details in the comments throughout this file. 
# If you are new to Quantarhei and/or Python and/or the simulated problem,
# it is recommended that you read the code, including the comments, as 
# a single text. All concepts of Quantarhei and
# the simulation model are explained when encountered, or a reference is
# provided. 
#
#
#
#
# Simulation Code
# ---------------
#
# First the components of the Quantarhei package used for the simulations 
# are imported to Python. Comments above the code lines explain what is 
# happening.
#

# Importing ...
# basic Quantarhei class to represent molecules
from quantarhei import Molecule

# basic Quantarhei class to represent 1D coordinate in molecules 
# (e.g. harmonic oscillator) 
from quantarhei import Mode

# basic Quantathei class to represent bath correlation function
from quantarhei import CorrelationFunction

# class representing the time interval on which we simulate
from quantarhei import TimeAxis

# 
# Comment on Quantarhei: internally, Quantarhei uses the units of cycles per
# femtosecond, or 2pi/fs, for energy. This is useful when simulating processes
# on femtosecond time scale. To define energies of optical transitions, one
# might find it useful (from a human perspective) to use something else. Below
# we use inverse centimeters (1/cm) which is a very common unit in
# spectroscopy. To tell Quantarhei that you enter values in other than the
# internal units, we use Python's context manager clause "with", and a special
# object "energy_units" below.
#
# context manager for physical units of energy is imported
from quantarhei import energy_units

from quantarhei import eigenbasis_of

# class representing Redfield relaxation rate matrix (and does its calculation)
from quantarhei.qm import RedfieldRateMatrix

from quantarhei import PopulationPropagator

# Here we import a very useful library for dealing with matrices
import numpy

# To produce figures, we import the pyplot from the matplotlib library
import matplotlib.pyplot as plt

from  quantarhei.core.managers import Manager

#
# First we create a molecule object representing the P state of the RC 
# and its charge transfer (CT) state. For the details of the mode, and its
# justification, see the article text (!!!!!)
#

print("\nQuantarhei version ",Manager().version)

# if "verbose" is set to True, the script prints some info on the screen
verbose = True

# if "show_plots" is set to True, the script will show plots
show_plots = False

if verbose:
    print("\nSingle photosynthetic reaction center fluorescence simulation")
    print(  "-------------------------------------------------------------\n")
    
#
# All definitions below are done with energy units of inverse 
# centimeters (1/cm). We use Python's context manager concept for this.
#
# Here we switch to the 1/cm units 
with energy_units("1/cm"):
    
    # Molecule named "RC-P" with three states (ground, CT and exciton state)
    # is created
    en = [0.0, 14500, 15500]
    m = Molecule("RC-P",en)    

    # Transition dipole moment from 0 -> 2 (exciton state) is set to 
    # a unit vector in the direction of x-axis (this choice is arbitrary)
    di = numpy.sqrt(23.0)*0.20819434
    m.set_dipole(0,2,[di,0.0,0.0])
    
    #
    # To represent CT coordinate, we add an explicite harmonic mode 
    # to the molecule.
    #
    
    # First the mode with frequency 10 1/cm is created
    qct = Mode(frequency=10.0)
    # and registered with the molecule
    m.add_Mode(qct)
    
    # In the CT state, the equilibrium coordinate is shifted towards
    # charge separation. In the ground state and the exciton, no shift
    # of the coordinate is assumed (see the text (!!!!) for justification.)
    qct.set_shift(1,1.0)
    
    # By default, the Mode we created is harmonic. The number of states
    # by which it is represented is set to 2 (which is a very small number).
    # Let us increase the number of states to represent the oscillators
    # in the three electronic states of the molecule to some chosen N[0], N[1]
    # and N[2] (see the text (!!!!) for the choice of values).
    #
    # we create an array of 3 zeros
    N = numpy.zeros(3,dtype=numpy.int)
    # set all of them to the value of 10
    N[0] = 10
    N[1] = 20
    N[2] = 20
    # loop from 0 to 2
    for i in range(3):
        # and set the number of vibrational states 
        qct.set_nmax(i,N[i])
    
    
    # Now we create a bath correlation function for the two electronic states.
    # This is done by creating a Python dictionary of parameters. These are:
    # - Type of the correlation function ftype = "OverdampedBrownian"
    # - Temperature T = 300 K
    # - Correlation time cortime = 100 fs
    # - Reorganization energy reorg = 100 1/cm (we are still in 
    #                                           the energy_units context)
    corfunc_ct_params = {"ftype":"OverdampedBrownian",
                         "T":300,"cortime":100,"reorg":100}
                         
    # We also need to specify a time interval on which the function is
    # defined. In principle, this will be the interval on which we will
    # calculate excited state dynamics. Quantarhei has an object which defines
    # the time interval. It is called TimeAxis. Arguments are 
    # 1) starting time, 2) Number of steps and 3) time step
    timeAxis = TimeAxis(0.0, 2000, 1.0)
    
    # Now we create a correlation function object
    corfunc_ct = CorrelationFunction(timeAxis,corfunc_ct_params)
    
    # The same is done below for a second correlation function, which represents
    # the interaction of the excitonic transition with the bath. Notice that
    # the reorganization energy is smaller here.
    corfunc_ex_params = {"ftype":"OverdampedBrownian",
                         "T":300,"cortime":100,"reorg":30}
    corfunc_ex = CorrelationFunction(timeAxis,corfunc_ex_params)
    
    # In addition to the above two system-bath interactions, we also have to
    # specify how the vibrational modes are damped in different electronic
    # states. We use the same system-bath interaction theory and require
    # correlation function to describe the coupling. Again we need two,
    # one for each excited electronic state. We do not need the relaxation
    # in the ground state for anything in our calculations, but because
    # the ground state is an isolated state, we can at least use it to
    # see how fast the relaxation is in a pure harmonic oscillator.
    corfunc_ct_q_params = {"ftype":"OverdampedBrownian",
                         "T":300,"cortime":100,"reorg":5.0}
    corfunc_ct_q = CorrelationFunction(timeAxis,corfunc_ct_q_params)
    
    corfunc_ex_q_params = {"ftype":"OverdampedBrownian",
                         "T":300,"cortime":100,"reorg":5.0}
    corfunc_ex_q = CorrelationFunction(timeAxis,corfunc_ex_q_params) 


    # Here we plot the bath correlation functions    
    if verbose and show_plots:
        print("Bath correlation functions:")
        # FIXME: give the plot a title and describe axes
        
        # CT correlation function
        # plot real part of the correlation function
        plt.plot(timeAxis.time,numpy.real(corfunc_ct.data),"-k")
        # plot imaginary part of the correlation function
        plt.plot(timeAxis.time,numpy.imag(corfunc_ct.data),"-r")
        #
        # Exciton correlation function
        # plot real part of the correlation function
        plt.plot(timeAxis.time,numpy.real(corfunc_ex.data),"-b")
        # plot imaginary part of the correlation function
        plt.plot(timeAxis.time,numpy.imag(corfunc_ex.data),"-g")
        plt.show()
    
    #    
    # Molecular transitions from the ground state will now be
    # assigned the bath correlation functions we defined above.
    # 
    # egcf = energy gap correlation function
    #
    # state 1 is the CT state
    m.set_egcf((0,1),corfunc_ct)        
    #
    # state 2 is the exciton state
    m.set_egcf((0,2),corfunc_ex)
    
    # 
    # Oscillator damping environments have to be specified in each electronic
    # state separately
    #
    m.set_mode_environment(mode=0,elstate=1,corfunc=corfunc_ct_q)
    m.set_mode_environment(mode=0,elstate=2,corfunc=corfunc_ex_q)
    
    # in the ground state we set the same environment as in the state 2
    ###m.set_mode_environment(mode=0,elstate=0,corfunc=corfunc_ex_q)
    

    # Now we will set adiabatic coupling between the CT state and the exciton
    # state 
    m.set_adiabatic_coupling(1,2,100.0)
    
    
#
# We will not set any more numbers, so we can leave the context manager
#    
    
    
# We ask the molecule object to create and return an object representing
# the Hamiltonian operator
hh = m.get_Hamiltonian()
#print(hh.data)

# Next we need transition dipole moment operator. It will be usefull 
# for setting initial condition of the simulation.
dd = m.get_TransitionDipoleMoment()
    
with energy_units("1/cm"):
#if True:
    # inspect the molecule by eyes
    if verbose:
        print("Molecule specification:")
        print(m)    

# The last ingredient missing before we can start calculating excited state
# dynamics is to determine ralaxation dynamics among the excited state
# levels. For this we need an object which holds information about the
# system-bath interaction in our system. Not suprisingly, the corresponding
# class is called SystemBathInteraction, and the Molecule class objects
# can return it upon request.
sbi = m.get_SystemBathInteraction(timeAxis)
    
# The sbi object now contains all information needed to calculate rates
# of energy transfer within the excited state manifold which is caused
# by the pure dephasing described by the correlation functions 
# "corfunc_ex" and "corfunc_ct". It also contains the damping of 
# vibrational modes which we prescribed by correlation functions
# "corfunc_ct_q" and "corfunc_ex_q" (!!!!)
    
# With the system-bath interaction at hand, we can calculate e.g. Redfield
# excitonic population relaxation matrix. 
# !!! This can take a while !!!
RRM = RedfieldRateMatrix(hh,sbi)
  
if verbose:
    print("Relaxation times")
    for k in range(hh.getDim(0)):
        for l in range(hh.getDim(0)):
            if (RRM.data[k,l] > 0) and ((l-k)==1):
                print(k,l,1.0/RRM.data[k,l]," fs")

#
# We are in position to calculate excited state dynamics of charge transfer
# process. 
#

# We set initial condition for the calculation. Sofar everything was
# specified in the site-basis. We will switch to exciton basis 


with eigenbasis_of(hh):
        
    # get the ground state thermal density matrix
    rho_eq = m.get_thermal_ReducedDensityMatrix()
    
    # we set initial condition by exciting the molecule with the polarization
    # along its transition dipole moment
    # 
    # # rho_ini = numpy.dot(dd.data[:,:,0],numpy.dot(rho_eq,dd.data[:,:,0]))
    rho_ini = rho_eq.excite_delta(dmoment=dd,epolarization=[1.0, 0.0, 0.0])
    

    # we normalize the population of excited state to 1.0
    rho_ini.normalize2(norm=1.0)
        
    # we need only population dynamics, so let us convert density matrix to
    # a vector of populations
    pop_ini = rho_ini.get_populations()


    # Now we set up a popagator for the population vector
    prop = PopulationPropagator(timeAxis,RRM)
        
    #
    # Here starts the intersting part! 
    #

    # We calculate population dynamics
    # !!! This can take a while !!!
    pop_in_time = prop.propagate(pop_ini)
    
    for n in range(pop_in_time.shape[1]):
        plt.plot(timeAxis.time,pop_in_time[:,n])
        
    lft_ct = m.get_electronic_natural_lifetime(1,epsilon_r=2.34)
    lft_ex = m.get_electronic_natural_lifetime(2,epsilon_r=2.34)
    
    print(lft_ct/1000000.0," ns")
    print(lft_ex/1000000.0," ns")
    
"""
    # We are after fluorescence. FluorescenceSpect class knows how to deal
    # with that. We submit the TimeAxis object and the Molecule object
    fl = FluorescenceSpect(timeAxis,m)
    
    #
    # calculate various spectra
    #
    
    # accumulated fluorescence calculated from equilibrium
    fspect = fl.accumulated_fluorescence()
    fspect.save_to_file("fluorecence.dat")
    
    # time correlated photon counting
    tcpc = fl.get_time_correlated_single_photon_counting(pop_in_time,
                                            spectral_range=(10000,18000),
                                            time_resolution=50000)
    tcpc.save_to_file("tcspc.dat")
    
"""

"""
# 
# Production calculations are based on scanning the position of CT state
# and then sampling the results randomly
#    
    
energy_of_CT_start = 9000
energy_of_CT_step = 100
energy_of_CT_Number_of_steps = 100

with energy_units("1/cm"):
    
    for n in range(energy_of_CT_Number_of_steps):
        energy_CT = energy_of_CT_start + energy_of_CT_step*n 
        m.elenergies[1] = energy_CT
        hh = m.get_Hamiltonian()
        
        #
        # Do the following objects change? 
        # 
        dd = m.get_TransitionDipoleMoment()
        sbi = m.get_SystemBathInteraction()
        
        rho_ini = numpy.dot(dd.data[:,:,0],numpy.dot(rho_eq,dd.data[:,:,0]))
        # we normalize the population of excited state to 1.0
        rho_ini = util.normalize_trace2(rho_ini,norm=1.0)
        pop_ini = rho_ini.get_population_vector()
        
        RRM = RedfieldRateMatrix(hh,sbi)
        prop = PopulationPropagator(timeAxis,RRM)
        pop_in_time = prop.propagate(pop_ini)
        
        # presumably we can use the same FluorescenceSpect object, because it
        # has an access
        # to the changed Molecule object
        fspect = fl.steady_state_fluorescence()
        tcpc = fl.get_time_correlated_photon_counting(pop_in_time,
                                            spectral_range=(10000,18000),
                                            time_resolution=50000)
                                            
        # append new results to the result file 
        
"""