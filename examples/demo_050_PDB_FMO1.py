# -*- coding: utf-8 -*-

_show_plots_ = False

print("""
***********************************************************      
*
*
*     PDB Extraction of FMO and Absorption Calculation
*
*
***********************************************************
""")


import numpy

from quantarhei.builders.pdb import PDBFile
from quantarhei.models.bacteriochlorophylls import BacterioChlorophyll

from quantarhei import Aggregate
from quantarhei import energy_units

import quantarhei as qr
    
#
# Read a PDB file
#
file = PDBFile("data_050_3eoj.pdb")
#file = PDBFile("data_050_3eni.pdb")
print("Loaded", file.linecount, "lines")

#
# Bacteriochlorophyll model will be used to extract molecules from PDB
#
bcl_model = BacterioChlorophyll(model_type="PDB")

#
# Exctract the molecules
#
molecules = file.get_Molecules(model=bcl_model)
names = []
for m in molecules:
    names.append(m.name)

names.sort()
print(names)

#
# Rename the molecules according to the standard from literature
# and exclude the BChl8
#
for_aggregate = []
naming_map = {"A371":"BChl1", "A372":"BChl2",
              "A373":"BChl3", "A374":"BChl4", "A375":"BChl5", 
              "A376":"BChl6", "A377":"BChl7", "A378":"BChl8"}

for name in names:
    for m in molecules:
        if m.name == name:
            m.set_name(naming_map[name])
            # all except for 378 go into aggregate
            if name != "378":
                for_aggregate.append(m)
             

        
#
# Create an new aggregate of the Bacteriochlorophylls without BChl8
#
agg = Aggregate(name="FMO", molecules=for_aggregate)


# Setting site energies according to literature
#
with energy_units("1/cm"):
    m = agg.get_Molecule_by_name("BChl1")
    m.set_energy(1, 12468.0)
    m = agg.get_Molecule_by_name("BChl2")
    m.set_energy(1, 12466.0)
    m = agg.get_Molecule_by_name("BChl3")
    m.set_energy(1, 12129.0)    
    m = agg.get_Molecule_by_name("BChl4")
    m.set_energy(1, 12410.0)
    m = agg.get_Molecule_by_name("BChl5")
    m.set_energy(1, 12320.0)
    m = agg.get_Molecule_by_name("BChl6")
    m.set_energy(1, 12593.0)
    m = agg.get_Molecule_by_name("BChl7")
    m.set_energy(1, 12353.0)
    
#
# Set resonance coupling by dipole-dipole method
#
agg.set_coupling_by_dipole_dipole(epsr=1.21)


# Setting spectral density to all molecules
tmax = qr.TimeAxis(0.0, 1000, 1.0)
param = dict(ftype="OverdampedBrownian", reorg=50,
            cortime=50, T=300)
with qr.energy_units("1/cm"):
    cfce = qr.CorrelationFunction(tmax, param)
    
for mol in agg.monomers:
    mol.set_transition_environment((0,1), cfce)


#
# Build the aggregate
#
agg.build()


#
# Now we can start simulations
#
H = agg.get_Hamiltonian()

with energy_units("1/cm"):
    print("\nExcited state Hamiltonian (energies in 1/cm):\n")
#    for i in range(1, H.dim):
#        for j in range(i,H.dim):
#            print(i,j,":",H.data[i,j])
    numpy.set_printoptions(precision=2, linewidth=100,
                           formatter={'all':lambda x: "%8.1f" % x})
    print(H.data[1:,1:])


agg.diagonalize()


calc = qr.AbsSpectrumCalculator(tmax, system=agg)
#with qr.energy_units("1/cm"):
rwa = agg.get_RWA_suggestion()
print("RWA frequency =", rwa)
calc.bootstrap(rwa=rwa)
abss = calc.calculate()
abss.normalize2()

if _show_plots_:
    with qr.energy_units("1/cm"):
        abss.plot(axis=[11500,13000,0,1.1])    


