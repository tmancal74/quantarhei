# -*- coding: utf-8 -*-

_show_plots_ = False
"""
    Demo of the Lindblad form
    
    First we create a Lindblad form for a purely electronic system,
    the we create electronic Lindblad form for a system with vibrational
    levels
    
"""
print("""
***********************************************************      
*
*
*          Lindblad form demo
*
*
***********************************************************
""")

from quantarhei import Manager
#
# FIXME: temporary fix for version 0.0.34
#
Manager().gen_conf.legacy_relaxation = True

#
# PURELY ELECTRONIC Aggregate of two molecules
#
from quantarhei import Molecule

m1 = Molecule([0.0, 1.0])
m2 = Molecule([0.0, 1.1])
m3 = Molecule([0.0, 1.2])

from quantarhei import Aggregate

agg = Aggregate([m1, m2, m3])
agg.build()

#
# Operator describing relaxation
#
from quantarhei.qm import Operator

HH = agg.get_Hamiltonian()
K = Operator(dim=HH.dim,real=True)
K.data[1,2] = 1.0

#
# System bath interaction with prescribed rate
#
from quantarhei.qm import SystemBathInteraction

sbi = SystemBathInteraction(sys_operators=[K], rates=(1.0/100,))
agg.set_SystemBathInteraction(sbi)

#
# Corresponding Lindblad form
#
from quantarhei.qm import LindbladForm

LF = LindbladForm(HH, sbi, as_operators=False)

print(LF.data[1,1,2,2])
print(LF.data[1,2,1,2])

#
# We can get it also from the aggregate
#


from quantarhei import TimeAxis

time = TimeAxis()

# time is not used here at all
LFa, ham = agg.get_RelaxationTensor(time, 
           relaxation_theory="electronic_Lindblad")
LFa.convert_2_tensor()

print(LFa.data[1,1,2,2])
print(LFa.data[1,2,1,2])


#
# VIBRONIC Aggregate of two molecules
#

m1v = Molecule([0.0, 1.0])
m2v = Molecule([0.0, 1.1])
m3v = Molecule([0.0, 1.2])

from quantarhei import Mode

mod1 = Mode(0.01)
mod2 = Mode(0.01)
mod3 = Mode(0.01)

Nvib = 2

m1v.add_Mode(mod1)
mod1.set_nmax(0, Nvib)
mod1.set_nmax(1, Nvib)
mod1.set_HR(1, 0.3)

m2v.add_Mode(mod2)
mod2.set_nmax(0, Nvib)
mod2.set_nmax(1, Nvib)
mod2.set_HR(1, 0.3)

m3v.add_Mode(mod3)
mod3.set_nmax(0, Nvib)
mod3.set_nmax(1, Nvib)
mod3.set_HR(1, 0.3)

vagg = Aggregate(molecules=[m1v, m2v, m3v])
vagg.build()

from quantarhei.qm import ProjectionOperator

K12 = ProjectionOperator(1, 2, dim=4)
K23 = ProjectionOperator(2, 3, dim=4)
ops = [K12, K23]
rates = [1.0/100.0, 1.0/150.0]

vsbi = SystemBathInteraction(sys_operators=ops, rates=rates)
vsbi.set_system(vagg)

from quantarhei.qm import ElectronicLindbladForm

ham = vagg.get_Hamiltonian()
print("Hamiltonian dimension: ", ham.dim)
print("Constructing electronic Lindblad form")
eLF = ElectronicLindbladForm(ham, vsbi, as_operators=True)
print("...done")

time = TimeAxis(0.0, 1000, 1.0)

from quantarhei.qm import ReducedDensityMatrixPropagator
vibprop = ReducedDensityMatrixPropagator(time, ham, eLF)

from quantarhei.qm import ReducedDensityMatrix
rho0 = ReducedDensityMatrix(dim=ham.dim)
Ni = vagg.vibindices[3][0]
rho0.data[Ni, Ni] = 1.0

print("Density matrix propagation: ")
rhot = vibprop.propagate(rho0)
print("...done")

print("\nThe same but NO vibrational modes")

m1v = Molecule([0.0, 1.0])
m2v = Molecule([0.0, 1.1])
m3v = Molecule([0.0, 1.2])

vagg2 = Aggregate(molecules=[m1v, m2v, m3v])
vagg2.build()

from quantarhei.qm import ProjectionOperator

K12 = ProjectionOperator(1, 2, dim=4)
K23 = ProjectionOperator(2, 3, dim=4)
ops = [K12, K23]
rates = [1.0/100.0, 1.0/150.0]

vsbi = SystemBathInteraction(sys_operators=ops, rates=rates)
vsbi.set_system(vagg2)

from quantarhei.qm import ElectronicLindbladForm

ham = vagg2.get_Hamiltonian()
print("Hamiltonian dimension: ", ham.dim)
print("Constructing electronic Lindblad form")
eLF = ElectronicLindbladForm(ham, vsbi, as_operators=False)
print("...done")

time = TimeAxis(0.0, 1000, 1.0)

from quantarhei.qm import ReducedDensityMatrixPropagator
vibprop = ReducedDensityMatrixPropagator(time, ham, eLF)

from quantarhei.qm import ReducedDensityMatrix
rho0 = ReducedDensityMatrix(dim=ham.dim)
Ni = 3
rho0.data[Ni, Ni] = 1.0

print("Density matrix propagation: ")
rhot2 = vibprop.propagate(rho0)
print("...done")

print("Trace over vibrations")
rhot1 = vagg.trace_over_vibrations(rhot) #, Nt=20)
print("...done")

#print(rhot2.data[Nt,:,:])
#print(rhot1.data)

if _show_plots_:
    rhot2.plot(show=False)
    rhot1.plot(how="--")

