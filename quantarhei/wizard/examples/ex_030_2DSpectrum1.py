# -*- coding: utf-8 -*-
"""

    Calculation of 2D spectra with effective lineshapes 


"""
import tempfile
import os

import numpy

import quantarhei as qr
import quantarhei.spectroscopy as spec

from quantarhei.qm.liouvillespace.evolutionsuperoperator import EvolutionSuperOperator


_show_plots_ = True

###############################################################################
#
# MODEL: Simple dimer of molecules
#
###############################################################################

with qr.energy_units("1/cm"):
    # two two-level molecules
    m1 = qr.Molecule([0.0, 12000.0])
    m2 = qr.Molecule([0.0, 12500.0])
    
    # transitions will have Gaussian lineshape with a width specified here
    m1.set_transition_width((0,1), 200.0)
    m2.set_transition_width((0,1), 200.0)

# we create an aggregate from the two molecules
agg = qr.Aggregate(molecules=[m1, m2])

# we set transition dipole moment orientations for the two molecules
m1.set_dipole(0,1,[1.0, 0.0, 0.0])
m2.set_dipole(0,1,[0.2, 0.9, 0.0])

# resonance coupling is set by hand
with qr.energy_units("1/cm"):
    agg.set_resonance_coupling(0,1, 100.0)

# we copy the aggregate before it is built. For the calculation of 2D 
# spectrum, we need to build the aggregate so that it contains two-exciton
# states. But those are irrelevant for single exciton excited state dynamics
# so we make two identical aggregates, one with single-excitons only, and
# one with two-excitons. 
#
# The copy is done with a workaround in which we
# save the to a temporary file and load it
with tempfile.TemporaryDirectory() as td:
    fname = os.path.join(td,"agg.qrp")
    qr.save_parcel(agg, fname)
    agg_2D = qr.load_parcel(fname)

# the aggregate is built with single exciton states only
agg.build(mult=1)

# we print its Hamiltonian to check everything is alright
H = agg.get_Hamiltonian()
with qr.energy_units("1/cm"):
    print(H)


###############################################################################
#
# EXCITED STATE DYNAMICS: Lindblad relaxation between eigenstates
#
###############################################################################

# time span of the excited state evolution (later t2 time of the 2D spectrum)
t2_axis = qr.TimeAxis(0.0, 100, 10.0)

K = qr.qm.ProjectionOperator(1,2,dim=H.dim)
rates = [1.0/200.0]
print(K)

sbi = qr.qm.SystemBathInteraction(sys_operators=[K], rates=rates)

L = qr.qm.LindbladForm(H, sbi)

print(L)

eUt = EvolutionSuperOperator(time=t2_axis, ham=H, relt=L)
eUt.set_dense_dt(10)
eUt.calculate()

with qr.eigenbasis_of(H):
    eUt.plot_element((2,2,2,2), show=False)
    eUt.plot_element((1,1,1,1), show=False)
    eUt.plot_element((1,1,2,2))
    eUt.plot_element((1,2,1,2))
    
    
###############################################################################
#
# 2D SPECTRUM: effective lineshape 2D spectrum
#
###############################################################################

# time axes of the propagation in t1 and t3 times

t1_axis = qr.TimeAxis(0.0, 100, 10.0)
t3_axis = qr.TimeAxis(0.0, 100, 10.0)

from quantarhei.spectroscopy.mocktwodcalculator \
    import MockTwoDSpectrumCalculator as TwoDSpectrumCalculator
    
from quantarhei.spectroscopy import X

calc = TwoDSpectrumCalculator(t1_axis, t2_axis, t3_axis)
print(H.rwa_energies)
calc.bootstrap(rwa=qr.convert(12200.0,"1/cm","int"))

agg_2D.build(mult=2)
agg_2D.diagonalize()

# laboratory settings
lab = qr.LabSetup()
lab.set_polarizations(pulse_polarizations=(X,X,X), detection_polarization=X)

T2 = 0.0
Uin = eUt.at(T2)

rho0 = agg_2D.get_DensityMatrix(condition_type="thermal",
                                      temperature=0.0)

# get Liouville pathways
pws = agg_2D.liouville_pathways_3T(ptype="R2g",
                                   eUt=Uin, ham=H, t2=T2, lab=lab)

pws = spec.order_by_amplitude(pws)
with qr.energy_units("1/cm"):
    print(len(pws))
    for pw in pws:
        print(pw)

calc.set_pathways(pws)

twod1 = calc.calculate_next()
twod1.set_t2(0.0)

print(numpy.max(twod1.data))
with qr.energy_units("1/cm"):
    twod1.plot(Npos_contours=10,              
              stype="total", spart="real")

