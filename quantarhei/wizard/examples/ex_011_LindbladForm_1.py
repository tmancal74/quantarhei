# -*- coding: utf-8 -*-

_show_plots_ = False

"""
    Demo of the Lindblad form
    
    First we create a Lindblad form for a purely electronic system,
    the we create electronic Lindblad form for a system with vibrational
    levels
    
"""
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
from quantarhei import Molecule

m1v = Molecule([0.0, 1.0])
m2v = Molecule([0.0, 1.1])
m3v = Molecule([0.0, 1.2])

from quantarhei import Mode

