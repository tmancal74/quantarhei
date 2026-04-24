
# coding: utf-8

# In[26]:
import os
import time

import quantarhei as qr

import matplotlib.pyplot as plt
plt.switch_backend('agg')

print(qr.Manager().version)


# In[27]:

pre_out = "out"
#frac = qr.load("fraction_40_4_CT_unbuilt.hdf5")
#frac_electronic = qr.load("fraction_40_4_CT_unbuilt.hdf5")

frac = qr.load_parcel(os.path.join(pre_out, "fraction_eff_40_4_CT_unbuilt.qrp"))
frac_electronic = qr.load_parcel(os.path.join(pre_out, "fraction_40_4_CT_unbuilt.qrp"))


# In[28]:


#
# Get molecules
#
PM = frac.get_Molecule_by_name("PM")
PL = frac.get_Molecule_by_name("PL")
BL = frac.get_Molecule_by_name("BL")

include_second_B = False 
vib_at_second_B = False 

if include_second_B:
    BM = frac.get_Molecule_by_name("BM")


# In[29]:


#
# Here we add vibrational modes
#

omega_B = 572
omega_P = 572   #572
with qr.energy_units("1/cm"):
    mod_PM = qr.Mode(omega_P)
    mod_PL = qr.Mode(omega_P)
    mod_BM = qr.Mode(omega_B)
    mod_BL = qr.Mode(omega_B)

Nmax_e = 1 
Nmax_g = 1 

HRF = 0.0111
vib_P = False 
vib_B = True 

if vib_P:
    PM.add_Mode(mod_PM)
    mod_PM.set_nmax(0, Nmax_g)
    mod_PM.set_nmax(1, Nmax_e)
    mod_PM.set_HR(1, HRF)

    PL.add_Mode(mod_PL)
    mod_PL.set_nmax(0, Nmax_g)
    mod_PL.set_nmax(1, Nmax_e)
    mod_PL.set_HR(1, HRF)

if vib_B:
    BL.add_Mode(mod_BL)
    mod_BL.set_nmax(0, Nmax_g)
    mod_BL.set_nmax(1, Nmax_e)
    mod_BL.set_HR(1, HRF)

    if include_second_B and vib_at_second_B:
        BM.add_Mode(mod_BM)
        mod_BM.set_nmax(0, Nmax_g)
        mod_BM.set_nmax(1, Nmax_e)
        mod_BM.set_HR(1, HRF)




# In[30]:


print(BL)


# In[31]:


#frac.save("fraction_45_2_vibrations_CT_unbuilt.hdf5")
qr.save_parcel(frac, os.path.join(pre_out, "fraction_45_2_vibrations_CT_unbuilt.qrp"))


# In[32]:


#
# Building the aggregate with approximations in vibrational manifold
#
frac.build(vibgen_approx="NPA", Nvib=2, mult=2) # NPA with Nvib=2 is equivalent to vibgen_approx="TPA"
frac_electronic.build(mult=2)


# In[33]:


print("Total number of states:", frac.Ntot)
print("Dim :", frac.get_Hamiltonian().data.shape)
print("Total number of electronic states:", frac.Nel, frac_electronic.Ntot)


# In[34]:


#
# Internals of the aggregates have to be converted to excitonic basis
#
import time
t1 = time.time()
frac.diagonalize()
t2 = time.time()
print("Diagonalized in ", t2-t1, "sec")
t3 = time.time()
frac_electronic.diagonalize()
print("Diagonalized in ", t3-t2, "sec")


# In[35]:


#
# Now we can look at the excitonic composition of all states
#
frac_electronic.report_on_expansion(1, N=4)

#
# State 1 is a clear P(-)
#


# In[36]:


frac_electronic.report_on_expansion(2, N=4)

#
# State 2 is a clear symmetric mixture of P(+) and B
#


# In[37]:


frac_electronic.report_on_expansion(3, N=4)

#
# State 3 is a clear anti-symmetric mixture of P(+) and B
#


# In[38]:


N1 = frac_electronic.Nbe[1]
print("Energies in 1/cm:")
print([qr.convert(frac_electronic.HH[i,i],"int","1/cm") for i in range(1, N1+1)])


# In[39]:


#
# Transition dipoles square
#
print("Transition dipoles square:")
print(frac_electronic.D2[1:N1+1,0])


# In[40]:


#frac.save("fraction_45_2_vibrations_CT_built.hdf5")
qr.save_parcel(frac, os.path.join(pre_out, "fraction_45_2_vibrations_CT_built.qrp"))


# In[41]:


frac.Nb


# In[42]:


frac.report_on_expansion(20)


# In[43]:


def criterium(lst):

    if lst[1] > 0.01:
        return True

    return False

with qr.energy_units("1/cm"):
    frac.exciton_report(start=6, criterium=criterium)


# # In[44]:
#
#
# #frac_B = qr.load("fraction_40_4_B_unbuilt.hdf5")
# #print(frac_B)
# frac_B = qr.load_parcel(os.path.join(pre_out, "fraction_40_4_B_unbuilt.qrp"))
#
#
# # In[45]:
#
#
# BL = frac_B.get_Molecule_by_name("BL")
# BM = frac_B.get_Molecule_by_name("BM")
#
#
# # In[46]:
#
# 
# omega_B = 572
#
# with qr.energy_units("1/cm"):
#     mod_BM = qr.Mode(omega_B)
#     mod_BL = qr.Mode(omega_B)
#
# Nmax_e = 2
# Nmax_g = 2
#
# HRF = 0.01
#
# vib_B = True
# if vib_B:
#     BL.add_Mode(mod_BL)
#     mod_BL.set_nmax(0, Nmax_g)
#     mod_BL.set_nmax(1, Nmax_e)
#     mod_BL.set_HR(1, HRF)
#
#     BM.add_Mode(mod_BM)
#     mod_BM.set_nmax(0, Nmax_g)
#     mod_BM.set_nmax(1, Nmax_e)
#     mod_BM.set_HR(1, HRF)
#
#
# # In[47]:
#
#
# #frac_B.save("fraction_45_2_vibrations_B_unbuilt.hdf5")
# qr.save_parcel(frac_B, os.path.join(pre_out, "fraction_45_2_vibrations_B_unbuilt.qrp"))
#
#
# # In[48]:
#
#
# #
# # Building the aggregate with approximations in vibrational manifold
# #
# frac_B.build(vibgen_approx="NPA", Nvib=2, mult=2) # NPA with Nvib=2 is equivalent to vibgen_approx="TPA"
#
#
# # In[49]:
#
#
# #
# # Internals of the aggregates have to be converted to excitonic basis
# #
# import time
# t1 = time.time()
# frac_B.diagonalize()
# t2 = time.time()
# print("Diaginalized in ", t2-t1, "sec")
# t3 = time.time()
#
#
# # In[50]:
#
#
# #frac_B.save("fraction_45_2_vibrations_B_built.hdf5")
# qr.save_parcel(frac_B, os.path.join(pre_out, "fraction_45_2_vibrations_B_built.qrp"))
# #qr.check_parcel("fraction_45_2_vibrations_B_built.qrp")
