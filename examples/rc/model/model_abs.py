# coding: utf-8
#
# Purely electronic model of the Reaction Center
#
# Calculations of absorption spectra with a realistic lineshape theory
# and with effective Gaussian lineshapes
#
#
#
#
# In[1]:

import os
import numpy

import quantarhei as qr
print(qr.Manager().version)


import matplotlib.pyplot as plt
plt.switch_backend('agg')

# In[2]:

pre_in = "in"
pre_out = "out"

# check if pre_out exists and is a directory
if not os.path.isdir(pre_out):
    try:
        os.makedirs(pre_out, exist_ok=True)
    except:
        raise Exception("Output directory name '"
                        +pre_out+"' does not represent a valid directory")


#
# Model from Jordanides at al. Ref. 1 is adjusted and extended by two CT states
#
#
jordanides = False 

if jordanides:
    offset = 0.0
    offset_P = 0.0 #485.0
    offset_P_M = offset_P + 0.0
    h_shift = 0.0
    sc_H = 1.0
    sc_P = 1.0
else:
    offset = 275
    offset_P = 400 #485.0
    offset_P_M = offset_P + 100.0
    h_shift = 85.0
    sc_H = 0.79
    sc_P = 0.75


#
# Molecules
#
with qr.energy_units("1/cm"):
    PM = qr.Molecule([0.0, 11610.0+offset_P_M], name="PM")
    PL = qr.Molecule([0.0, 11610.0+offset_P], name="PL")
    BM = qr.Molecule([0.0, 12220.0+offset], name="BM")
    BL = qr.Molecule([0.0, 12370.0+offset], name="BL")
    HL = qr.Molecule([0.0, 13020.0+offset-h_shift], name="HL")
    HM = qr.Molecule([0.0, 13150.0+offset+h_shift], name="HM")

    # CT states are effectively represented as "new molecules" in the system
    PCT_M = qr.Molecule([0.0, 15200], name="PCT1")
    PCT_L = qr.Molecule([0.0, 13550], name="PCT2")  # 13500

#
# Transition dipole moment from Ref. 1 are scaled
#
dPM = numpy.array([ 0.8546,  0.5051,  0.1206])*sc_P
dPL = numpy.array([-0.9649, -0.0250,  0.2613])*sc_P
dHM = numpy.array([ 0.2749, -0.3694, -0.8877])*sc_H
dHL = numpy.array([ 0.0452, -0.9672, -0.2498])*sc_H


PM.set_dipole(0,1, dPM)
PL.set_dipole(0,1, dPL)
BL.set_dipole(0,1, [ 0.7782,  0.5332,  0.3317])
BM.set_dipole(0,1, [-0.9681,  0.1107,  0.2249])
HL.set_dipole(0,1, dHL)
HM.set_dipole(0,1, dHM)

#
# CT states are dark
#
PCT_M.set_dipole(1, 0, [0.0, 0.0, 0.0])
PCT_L.set_dipole(1, 0, [0.0, 0.0, 0.0])

molecules = [PM, PL, BM, BL, HL, HM, PCT_M, PCT_L]

# saving molecules without environment
qr.save_parcel(molecules, os.path.join(pre_out,"molecules.qrp"))

#
# Here we build the RC as an aggregate of molecules
#
mol3 = [PM, PL, BM]
agg = qr.Aggregate(molecules=mol3)


#
# Exciton interaction matrix
#

# values from Ref. 1
JP_77K_Jordanides = 575.0
JP_77K = JP_77K_Jordanides

#
# Fitted values of the model with CT states
# starting values of the manual search of best parameters are
# taken from Ref. 2
#
if jordanides:
    JP = 395 #JP_77K
    XCT_M = 0.0
    XCT_L = 0.0
    YCT = 0.0
else:
    JP = 690 #575
    XCT_M = 905 #1400
    XCT_L = 755
    YCT = 550 #350

# Factor of three is just to experiment with
PB_1 = -104.0
PB_2 = -94.0

LCT = 0
MCT = 0

# the interaction matrix is taken from
J_Matrix = numpy.array([
    [   0.0,    JP, -16.0, PB_1,  19.9, -4.8, XCT_M, YCT],
    [    JP,   0.0, PB_2,    2.8,  -6.8, 18.0, YCT, XCT_L],
    [ -16.0, PB_2,   0.0,   19.3,  -7.5, 95.8, MCT, LCT],
    [ PB_1,   2.8,  19.3,    0.0, 123.1, -7.9, LCT, MCT],
    [  19.9,  -6.8,  -7.5,  123.1,   0.0,  3.9, 0.0, 0.0],
    [  -4.8,  18.0,  95.8,   -7.9,   3.9,  0.0, 0.0, 0.0],
    [   XCT_M,   YCT,   MCT,    LCT,   0.0,  0.0, 0.0, 0.0],
    [   YCT,   XCT_L,   LCT,    MCT,   0.0,  0.0, 0.0, 0.0]
])

with qr.energy_units("1/cm"):
    agg.set_resonance_coupling_matrix(J_Matrix[0:3,0:3])

#agg.save("RC_Model_40_4_adjusted_CT_no_environment_unbuilt.hdf5")
qr.save_parcel(agg, os.path.join(pre_out,
                                 "RC_Model_40_4_adjusted_CT_no_environment_unbuilt.qrp"))


# In[3]:


# check that units were set correctly
rc = agg.resonance_coupling[1,0]
with qr.energy_units("1/cm"):
    print(qr.convert(rc, "int"))

with qr.energy_units("1/cm"):
    print(agg.get_resonance_coupling(1,0))


# In[4]:


# Bath correlation function
time = qr.TimeAxis(0.0, 1000, 1.0)

cfA_params = dict(ftype="OverdampedBrownian",
                  reorg=190, cortime=80, T=77, matsubara=100)
cfH_params = dict(ftype="OverdampedBrownian",
                  reorg=200, cortime=100, T=77, matsubara=100)
cfP_params = dict(ftype="OverdampedBrownian",
                  reorg=700, cortime=120, T=77, matsubara=100)
cfCT_params = dict(ftype="OverdampedBrownian",
                  reorg=3600, cortime=20, T=77, matsubara=200)


with qr.energy_units("1/cm"):
    cfA = qr.CorrelationFunction(time, cfA_params)
    cfH = qr.CorrelationFunction(time, cfH_params)
    cfP = qr.CorrelationFunction(time, cfP_params)
    cfCT = qr.CorrelationFunction(time, cfCT_params)

PM.set_transition_environment((0,1), cfP)
PL.set_transition_environment((0,1), cfP)
BM.set_transition_environment((0,1), cfA)
BL.set_transition_environment((0,1), cfA)
HL.set_transition_environment((0,1), cfH)
HM.set_transition_environment((0,1), cfH)
PCT_M.set_transition_environment((0,1), cfCT)
PCT_L.set_transition_environment((0,1), cfCT)

agg.build(mult=2)
#agg.save("RC_Model_40_4_adjusted_CT_no_vibrations_built.hdf5")
qr.save_parcel(agg, os.path.join(pre_out,
    "RC_Model_40_4_adjusted_CT_no_vibrations_built.qrp"))

# In[5]:

#
#  Refitted model of the Reaction Center using effective Gaussian lineshapes
#
#
#
#
# "Environment" modelled by dressing the states
#
molecules_eff = qr.load_parcel(os.path.join(pre_out,"molecules.qrp"))

agg_eff = qr.Aggregate(molecules=molecules_eff)
with qr.energy_units("1/cm"):
    agg_eff.set_resonance_coupling_matrix(J_Matrix)

PMe = molecules_eff[0]
PLe = molecules_eff[1]
BMe = molecules_eff[2]
BLe = molecules_eff[3]
HMe = molecules_eff[4]
HLe = molecules_eff[5]
PCT_Me = molecules_eff[6]
PCT_Le = molecules_eff[7]

with qr.energy_units("1/cm"):
    ee = PMe.get_energy(1)
    PMe.set_energy(1,ee-80.0)
    ee = PLe.get_energy(1)
    PLe.set_energy(1,ee-80.0)
    ee = BMe.get_energy(1)
    BMe.set_energy(1,ee-85.0)
    ee = BLe.get_energy(1)
    BLe.set_energy(1,ee-85.0)
    ee = HMe.get_energy(1)
    HMe.set_energy(1,ee-75.0)
    ee = HLe.get_energy(1)
    HLe.set_energy(1,ee-75.0)

    ee = PCT_Me.get_energy(1)
    PCT_Me.set_energy(1,ee+230)
    ee = PCT_Le.get_energy(1)
    PCT_Le.set_energy(1,ee+230)

PMe.set_transition_width((0,1), qr.convert(630,"1/cm", "int"))
PLe.set_transition_width((0,1), qr.convert(630,"1/cm", "int"))
BMe.set_transition_width((0,1), qr.convert(180,"1/cm", "int"))
BLe.set_transition_width((0,1), qr.convert(180,"1/cm", "int"))
HMe.set_transition_width((0,1), qr.convert(155,"1/cm", "int"))
HLe.set_transition_width((0,1), qr.convert(155,"1/cm", "int"))
PCT_Me.set_transition_width((0,1), qr.convert(800,"1/cm", "int"))
PCT_Le.set_transition_width((0,1), qr.convert(800,"1/cm", "int"))

dPMe = numpy.array([ 0.8546,  0.5051,  0.1206])*0.76
dPLe = numpy.array([-0.9649, -0.0250,  0.2613])*0.76
dHMe = numpy.array([ 0.2749, -0.3694, -0.8877])*0.68
dHLe = numpy.array([ 0.0452, -0.9672, -0.2498])*0.68

PMe.set_dipole(0,1,dPMe)
PLe.set_dipole(0,1,dPLe)
HMe.set_dipole(0,1,dHMe)
HLe.set_dipole(0,1,dHLe)

# we save the effective model
qr.save_parcel(agg_eff, os.path.join(pre_out,
    "RC_eff_Model_40_4_adjusted_CT_no_environment_unbuilt.qrp"))

agg_eff.build(mult=1)


# In[6]:




#RT = agg.get_RelaxationTensor(time, relaxation_theory="standard_Redfield", secular_relaxation=True)
rrm = agg.get_RedfieldRateMatrix()


# In[7]:


print("Relaxation time (2 -> 1) :", 1.0/rrm.data[1,2])
print("Relaxation time (3 -> 2) :", 1.0/rrm.data[2,3])
print("Relaxation time (3 -> 1) :", 1.0/rrm.data[1,3])

# TEST (put here the energies and temperature)
E2 = 2.24165620051
E1 = 2.13494501445
kbT = 0.01008086552556262
print("Relaxation time ratio   :", rrm.data[2,1]/rrm.data[1,2])
print("... to be compared with :", numpy.exp(-(E2-E1)/kbT))


# In[8]:


rwa = agg.get_RWA_suggestion()
with qr.energy_units("1/cm"):
    print(qr.convert(rwa,"int"))


# In[9]:

# absorption from effective theory
from quantarhei import LabSetup
from quantarhei.utils.vectors import X, Y, Z

lab = LabSetup()
lab.set_polarizations(pulse_polarizations=[X,X,X], detection_polarization=X)

agg_eff.diagonalize()

print("\nEffetive model exciation energies:")
print("Energies in 1/cm:")
N1 = agg_eff.nmono
print([qr.convert(agg_eff.HH[i,i],"int","1/cm") for i in range(1, N1+1)])
print("")

mabsc = qr.MockAbsSpectrumCalculator(time, system=agg_eff)
rho0 = agg_eff.get_DensityMatrix(condition_type="thermal", temperature=0.0)
ham = agg_eff.get_Hamiltonian()

pthways = agg_eff.liouville_pathways_1(lab=lab, ham=ham, etol=1.0e-5,
                                       verbose=0)

mabsc.bootstrap(rwa=qr.convert(10000.0,"1/cm","int"),
              shape="Gaussian")

mabsc.set_pathways(pthways)

abs1 = mabsc.calculate(raw=False)
abs1.normalize2()

absc = qr.AbsSpectrumCalculator(time, system=agg)


# In[10]:


absc.bootstrap(rwa)


# In[11]:


abss = absc.calculate()
#absexp = qr.load("bas_77K.hdf5") #("DATA/bas_77K.hdf5")
absexp = qr.load_parcel(os.path.join(pre_in, "bas_77K.qrp"))

absexp.normalize()
absexp.subtract(0.086)
absexp.normalize()
abss.normalize2() #norm=0.53)


# In[29]:


#with qr.energy_units("nm"):
#    abss.plot(axis=[650, 1000, 0, 0.7], show=False)
#    absexp.plot()

plt.figure(0)
with qr.energy_units("1/cm"):
    #abss.plot(axis=[10500, 15000, 0, 1.1], show=False)
    abs1.plot(axis=[10500, 15000, 0, 1.1], show=False)
    absexp.plot(show=True)
    absexp.savefig(os.path.join(pre_out, "abs_full.png"))

# in a Notebook, it seems that the figure shows itself always when we leave the cell


# In[30]:


N1 = agg.Nbe[1]
print("Energies in 1/cm:")
print([qr.convert(agg.HH[i,i],"int","1/cm") for i in range(1, N1+1)])


# In[31]:


agg.diagonalize()


# In[32]:


# exciton report
agg.exciton_report(Nrep=8)


# In[33]:


agg.report_on_expansion(2)


# In[34]:


N1 = agg.Nbe[1]
print("Energies in 1/cm:")
print([qr.convert(agg.HH[i,i],"int","1/cm") for i in range(1, N1+1)])


# In[35]:


print("Transition dipoles square:")
print(agg.D2[1:N1+1,0])

#
# ## Fractional model
#
# Remove both H and BL

# In[36]:


#
# Get components of the fractional model
#
indices_of_components = []
names_of_components = ["PM", "PL", "BM"] # , "BL","PCT1", "PCT2"] # "HL", "HM", "BL", , "BCT"
components = []
for name in names_of_components:
    indx = agg.get_Molecule_index(name)
    mol = agg.get_Molecule_by_name(name)
    #if name == "BM":
    #    mol.elenergies[1] = mol.elenergies[1] + 0.1
    indices_of_components.append(indx)
    components.append(mol)
print("Indices of selected molecules: ", indices_of_components)


# In[37]:


#
# Coupling matrix
#
Ni = len(indices_of_components)
Jfm = numpy.zeros((Ni, Ni), dtype=qr.REAL)
k_1 = 0
for i_1 in indices_of_components:
    k_2 = 0
    for i_2 in indices_of_components:
        Jfm[k_1, k_2] = agg.resonance_coupling[i_1, i_2]
        k_2 += 1
    k_1 += 1


# In[38]:


#
# Fractional aggregate
#
frac = qr.Aggregate(components)
frac.set_resonance_coupling_matrix(Jfm)


# In[39]:


fix_dipole = False
if fix_dipole:
    BM_fix_dipole = frac.get_Molecule_by_name("BM")
    dip = BM_fix_dipole.get_dipole(0, 1)
    nrm = qr.norm(dip)
    dip2 = qr.normalize2(dip, norm=numpy.sqrt(2.0)*nrm)
    BM_fix_dipole.set_dipole(0, 1, dip2)


# In[40]:


#frac.save("fraction_40_4_CT_unbuilt.hdf5")
qr.save_parcel(frac, os.path.join(pre_out,"fraction_40_4_CT_unbuilt.qrp"))


# In[41]:


frac.build()


# In[42]:





absc2 = qr.AbsSpectrumCalculator(time, system=frac)
absc2.bootstrap(rwa)
abss2 = absc2.calculate()
#absexp2 = qr.load("bas_77K.hdf5")
absexp2 = qr.load_parcel(os.path.join(pre_in, "bas_77K.qrp"))

absexp2.normalize()
absexp2.subtract(0.086)
absexp2.normalize()
abss2.normalize2() #norm=0.53)

plt.figure(1)
with qr.energy_units("1/cm"):
    abss2.plot(axis=[10500, 15000, 0, 1.1], show=False)
    absexp2.plot(show=True)
    absexp2.savefig(os.path.join(pre_out, "abs_frac.png"))


# In[43]:





frac.diagonalize()


# In[44]:


frac.report_on_expansion(3)


# In[45]:


HH = frac.get_Hamiltonian()
with qr.eigenbasis_of(HH):
    with qr.energy_units("1/cm"):
        print([HH.data[i,i] for i in range(1,frac.nmono)])


# In[46]:


#
# Get components of the fractional model
#
indices_of_components = []
names_of_components = ["PM", "PL", "BM", "BL","PCT1", "PCT2"] #["BM", "BL"] # "HL", "HM", "BL", , "BCT"
names_of_components3 = ["PM", "PL", "BL"] 
components = []
for name in names_of_components3:
    indx = agg_eff.get_Molecule_index(name)
    mol = agg_eff.get_Molecule_by_name(name)
    #if name == "BM":
    #    mol.elenergies[1] = mol.elenergies[1] + 0.1
    indices_of_components.append(indx)
    components.append(mol)
print("Indices of selected molecules: ", indices_of_components)

# In[47]:


#
# Fractional aggregate
#
frac_eff = qr.Aggregate(components)
frac_eff.set_resonance_coupling_matrix(Jfm)


# In[48]:


#frac_B.save("fraction_40_4_B_unbuilt.hdf5")
qr.save_parcel(frac_eff, os.path.join(pre_out,
                         "fraction_eff_40_4_CT_unbuilt.qrp"))

frac_eff.build()
frac_eff.diagonalize()

mabsc2 = qr.MockAbsSpectrumCalculator(time, system=frac_eff)
rho0 = frac_eff.get_DensityMatrix(condition_type="thermal", temperature=0.0)
ham = frac_eff.get_Hamiltonian()

pthways = frac_eff.liouville_pathways_1(lab=lab, ham=ham, etol=1.0e-5,
                                       verbose=0)

mabsc2.bootstrap(rwa=qr.convert(10000.0,"1/cm","int"),
              shape="Gaussian")

mabsc2.set_pathways(pthways)

abs2 = mabsc2.calculate(raw=False)
abs2.normalize2()

plt.figure(2)
with qr.energy_units("1/cm"):
    #abss2.plot(axis=[10500, 15000, 0, 1.1], show=False)
    abs2.plot(axis=[10500, 15000, 0, 1.1], show=False)
    absexp2.plot(show=False)
    absexp2.savefig(os.path.join(pre_out, "abs_frac_eff.png"))
