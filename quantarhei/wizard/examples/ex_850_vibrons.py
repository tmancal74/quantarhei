# -*- coding: utf-8 -*-
import numpy

import quantarhei as qr
from quantarhei import printlog as print
from quantarhei.utils.vectors import X

_show_plots_ = False

# energy of molecule one
E1 = 12500.0
# energy gap to molecule two
Edelta = 100.0
# coupling between the two molecules
JJ = 30.0

# frequency of the vibrational mode
omega = 110.0
# Huan-Rhys factor
HR = 0.01

# transition width
width = 80

E2 = E1 + Edelta
print("")
print("Molecular dimer")
print("E1:", E1,"1/cm")
print("E2:", E2,"1/cm (delta =",Edelta,")")

with qr.energy_units("1/cm"):
    mol1 = qr.Molecule([0.0, E1])
    mol1.set_dipole(0,1,[1.0, 0.0, 0.0])
    mol1.set_transition_width((0,1), width)
    
    mod = qr.Mode(omega)
    mol1.add_Mode(mod)
    mod.set_nmax(0, 2)
    mod.set_nmax(1, 2)
    mod.set_HR(1, HR)
    
    mol2 = qr.Molecule([0.0, E2])
    mol2.set_dipole(0,1,numpy.array([1.0, 1.0, 0.0])/numpy.sqrt(2.0))
    mol2.set_transition_width((0,1), width)
    
    
    agg = qr.Aggregate(molecules=[mol1, mol2])
    agg.set_resonance_coupling(0,1,JJ)

agg.build(mult=2)
agg.diagonalize()

H = agg.get_Hamiltonian()

#with qr.energy_units("1/cm"):
#    print(H)

lab = qr.LabSetup()
lab.set_polarizations(pulse_polarizations=[X,X,X],
                      detection_polarization=X)

#print(qr.convert(agg.Wd,"int","1/cm"))

t1 = qr.TimeAxis(0.0, 1000, 5.0)
t3 = qr.TimeAxis(0.0, 1000, 5.0)
t2 = qr.TimeAxis(0.0, 2, 10.0)
calc = qr.MockTwoDResponseCalculator(t1,t2,t3)
with qr.energy_units("1/cm"):
    calc.bootstrap(rwa=12500.0)

eUt = qr.qm.EvolutionSuperOperator(t2, H)
eUt.set_dense_dt(20)

#
# We calculate evolution superoperator
#
eUt.calculate(show_progress=False)


olow_cm = omega-30.0
ohigh_cm = omega+30.0
olow = qr.convert(olow_cm, "1/cm", "int")
ohigh = qr.convert(ohigh_cm, "1/cm", "int")

sel_1 = [["omega2",[olow, ohigh]]]
sel_2 = [["omega2",[-ohigh, -olow]]]
#sel_1 = None
#sel_2 = None

pways = dict()
resp_plus = calc.calculate_one_system(t2.data[0], agg, eUt, lab, pways=pways,
                                      selection=sel_1,
                                      dtol=0.0001)

resp_mins = calc.calculate_one_system(t2.data[0], agg, eUt, lab, pways=pways,
                                      selection=sel_2,
                                      dtol=0.0001)

with qr.energy_units("1/cm"):
    twod = resp_plus.get_TwoDSpectrum()
    if _show_plots_:
        twod.plot(window=[12000,13000,12000,13000])

    twod = resp_mins.get_TwoDSpectrum()    
    if _show_plots_:
        twod.plot(window=[12000,13000,12000,13000])
    
    pw = []
    print(len(pways.keys()))
    for key in pways.keys():
        for p in pways[key]:
            pw.append(p)
            print(p)

len(pways[key])
print("%%%%%%%%%%%%%%%%")

p_way = pw[5]
calc.set_pathways([p_way])

resp = calc.calculate()
twod = resp.get_TwoDSpectrum()

if _show_plots_:
    with qr.energy_units("1/cm"):
        print(p_way)
        twod.plot(window=[12000,13000,12000,13000])

s = [2,3,1]
for ss in s:
    print("state",ss)
    agg.report_on_expansion(ss)

print(agg.DD[0,3])
print(agg.DD[0,2])
print(agg.DD[3,1])
print(agg.DD[2,1])