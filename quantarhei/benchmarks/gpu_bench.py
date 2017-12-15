# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 09:12:29 2017

@author: tomas
"""
import time
import os

import numpy

import quantarhei as qr

manager = qr.Manager()

#
# BENCHMARK SETTINGS
#


Nthreads = 4
use_mpi = False

fix_seed = True

N_molecules = 4 


#
# build an elementary model to propagate
#
qr.log_report("")
qr.log_report("Quantarhei benchmark calculation (task no. 001)",
              incr_indent=2)
qr.log_report("")

t_start = time.time()

benchmark_report = {}
benchmark_report["benchmark_id"] = "001"
benchmark_report["N_molecules"] = N_molecules


os.environ["OMP_NUM_THREADS"] = str(Nthreads)
benchmark_report["Nthreads"] = Nthreads


benchmark_report["gpu_acceleration"] = manager.num_conf.gpu_acceleration


if fix_seed:
    seed = 1234
    numpy.random.seed(seed)
    benchmark_report["seed"] = seed
else:
    benchmark_report["seed"] = False


benchmark_report["use_mpi"] = use_mpi
if not use_mpi:
    benchmark_report["mpi_comm_size"] = 1
else:
    benchmark_report["mpi_comm_size"] = -1

# Group of two-level molecules

ex_energies = numpy.zeros(N_molecules)
ex_energies[:] = 12000.0
mols = []
previous_position = numpy.zeros(3, dtype=qr.REAL)
dist_max = 10.0
dist_min = 5.0
placed_molecules = []
dipole_length = 6.0 

with qr.energy_units("1/cm"):
    qr.log_report("Generating molecules")
    qr.log_report("--------------------")
    for i_m in range(N_molecules):
        mol = qr.Molecule([0.0, ex_energies[i_m]])
        mols.append(mol)
        qr.log_info("Molecule:", i_m)
        
        # Randomly oriented transition dipole moment
        vec = numpy.random.rand(3) - 0.5*numpy.ones(3)
        vec = qr.normalize2(vec, dipole_length)
        mol.set_dipole(0, 1, vec)
        qr.log_info("Transition dipole moment:", vec)
        
        # Randomly placing the molecule
        placed = False
        while not placed:
            pos = numpy.random.rand(3) - 0.5*numpy.ones(3)
            pos = qr.normalize2(pos, dist_max)
            pos_n = numpy.array(pos) + previous_position
            # check distance to all molecules
            dist_OK = True
            for pmol in placed_molecules:
                rd = pos_n - pmol.position
                dist = numpy.sqrt(numpy.dot(rd, rd))
                if (dist < dist_min):
                    dist_OK = False
                    #print("Refusing distance: ", dist)
                else:
                    #print("New distance:", dist)
                    pass
            if dist_OK:
                mol.position = pos_n
                placed_molecules.append(mol)
                placed = True
                previous_position = pos_n
                qr.log_info("Placing molecule no.", i_m, "at:", mol.position)
              
        
    
# Vibrational modesif True:
Ne_max = 2
Ng_max = 2
vomega = 100
hr_fac = 0.01
with qr.energy_units("1/cm"):
    for i_m in range(N_molecules):
        mod = qr.Mode(vomega)
        mols[i_m].add_Mode(mod)
        mod.set_nmax(0, Ng_max)
        mod.set_nmax(1, Ne_max)
        mod.set_HR(1, hr_fac)
        
        
# Bath correlation functions
timea = qr.TimeAxis(0.0, 1000, 1.0)
cfpar = dict(ftype="OverdampedBrownian", 
             reorg=30, cortime=100, T=300, matsubara=100)
with qr.energy_units("1/cm"):
    cf = qr.CorrelationFunction(timea, cfpar)
    
# set bath correlation functions to the molecules
for i_m in range(N_molecules):
    mols[i_m].set_transition_environment((0, 1), cf)
    
# aggregate of molecules
agg = qr.Aggregate(mols)

agg.set_coupling_by_dipole_dipole()

# Building the aggregate
qr.log_report("Building aggregate")
agg.build()
qr.log_report("...done")

qr.log_detail("Resonance coupling matrix: ")
qr.log_detail(qr.convert(agg.resonance_coupling, "int", "1/cm"),
             use_indent=False)

# Dimension of the problem
HH = agg.get_Hamiltonian()
Nr = HH.dim
qr.log_detail("Hamiltonian has a rank:", Nr)

benchmark_report["Dimension"] = Nr
    
qr.log_report("Calculating Relaxation tensor:")
t1 = time.time()
(RT, ham) = agg.get_RelaxationTensor(timea, 
                                     relaxation_theory="standard_Redfield",
                                     as_operators=True)
t2 = time.time()
qr.log_report("...done in", t2-t1, "sec")

benchmark_report["standard_Redfield_timedependent_False"] = t2-t1

#prop = agg.get_ReducedDensityMatrixPropagator(timea, 
#                                       relaxation_theory="standard_Redfield",
#                                       as_operators=True )

prop = qr.ReducedDensityMatrixPropagator(timea, Ham=ham, RTensor=RT)


rho0 = agg.get_DensityMatrix(condition_type="impulsive_excitation", 
                             temperature=0.0)
rho0.normalize2()

qr.log_report("Propagating density matrix:")
t1 = time.time()
rhot = prop.propagate(rho0)
t2 = time.time()
qr.log_report("...done in", t2-t1, "sec")

benchmark_report["rdm_propagation_RTtimedependent_False"] = t2-t1

t_end = time.time()
qr.log_report()
qr.log_report("Benchmark calculation finished in ", t_end-t_start, "sec")
qr.log_report()

print(benchmark_report)


    
