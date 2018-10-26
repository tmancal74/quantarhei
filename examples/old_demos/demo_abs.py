# -*- coding: utf-8 -*-

import quantarhei as qr


with qr.energy_units("1/cm"):
    m1 = qr.Molecule([0.0, 10000.0])
    m1.set_dipole(0,1,[1.0, 0.0, 0.0])
    m2 = qr.Molecule([0.0, 11000.0])
    m2.set_dipole(0,1,[1.0, 0.0, 0.0])
    
agg = qr.Aggregate(molecules=[m1, m2])

time = qr.TimeAxis(0.0, 10000, 1.0)
cpar = dict(ftype="OverdampedBrownian", cortime=30, reorg=200, T=300)

with qr.energy_units("1/cm"):
    cf = qr.CorrelationFunction(time, cpar)

m1.set_transition_environment((0,1), cf)
m2.set_transition_environment((0,1), cf)

agg.build()

asc = qr.AbsSpectrumCalculator(time, system=agg)
asc.bootstrap(rwa=qr.convert(10050,"1/cm", "int"))

abs1 = asc.calculate(raw=False)
abs1.normalize2()

with qr.energy_units("1/cm"):
    abs1.plot(axis=[9000,12000,0.0,1.01], color="b", show=False)
    


with qr.energy_units("1/cm"):
    m3 = qr.Molecule([0.0, 9900.0])
    m3.set_dipole(0,1,[1.0, 0.0, 0.0])
    m4 = qr.Molecule([0.0, 10900.0])
    m4.set_dipole(0,1,[1.0, 0.0, 0.0])
    
    
#agg.clean()
m3.set_transition_width((0,1), qr.convert(380,"1/cm", "int"))
m4.set_transition_width((0,1), qr.convert(380,"1/cm", "int"))

agg2 = qr.Aggregate(molecules=[m3, m4])
agg2.build()

from quantarhei import LabSetup
from quantarhei.utils.vectors import X, Y, Z

lab = LabSetup()
lab.set_polarizations(pulse_polarizations=[X,X,X], detection_polarization=X)


agg2.diagonalize()

#
#  Absorption Spectrum by pathway method
#

mac = qr.MockAbsSpectrumCalculator(time, system=agg2)

rho0 = agg2.get_DensityMatrix(condition_type="thermal", temperature=0.0)
ham = agg2.get_Hamiltonian()

pthways = agg2.liouville_pathways_1(lab=lab, ham=ham, etol=1.0e-5,
                                       verbose=0) 

mac.bootstrap(rwa=qr.convert(10000.0,"1/cm","int"), 
              shape="Gaussian")

mac.set_pathways(pthways)

abs1 = mac.calculate(raw=False)
abs1.normalize2()

with qr.energy_units("1/cm"):
    abs1.plot(axis=[9000,12000,0.0,1.01], color="r")