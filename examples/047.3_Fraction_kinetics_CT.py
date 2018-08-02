
# coding: utf-8

# In[1]:


import quantarhei as qr
import quantarhei.spectroscopy as spec
import numpy
import time

print(qr.Manager().version)

import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


#
# This is needed in version 0.0.36 for propagation with Lindblad form
#
qr.Manager().gen_conf.legacy_relaxation = True


# In[3]:


agg = qr.load("fraction_45_2_vibrations_CT_unbuilt.hdf5")
agg2 = qr.load("fraction_45_2_vibrations_CT_unbuilt.hdf5")
agg_el = qr.load("fraction_40_4_CT_unbuilt.hdf5")


# In[4]:


#
# Aggregate with vibrational states, Hamiltonian generated up to single exciton states in a Two-particle approximation
# This object is used for calculation of dynamics
#
agg.build(mult=1, vibgen_approx="TPA")

#
# width of the spectral features in 1/cm
#
wincm = 500

#
# agg2 object will be used for calculation of 2D spectra
# so it needs to know how to represent the spectra (what
# width it should have), and it has to be built with
# two-exciton states (in two-particle approximation)
#
width = qr.convert(wincm, "1/cm", "int")
PM = agg2.get_Molecule_by_name("PM")
PM.set_transition_width((0,1), width)
PL = agg2.get_Molecule_by_name("PL")
PL.set_transition_width((0,1), width)

width = qr.convert(wincm, "1/cm", "int")
BM = agg2.get_Molecule_by_name("BM")
BM.set_transition_width((0,1), width)
BL = agg2.get_Molecule_by_name("BL")
BL.set_transition_width((0,1), width)

width = qr.convert(wincm, "1/cm", "int")
PCT1 = agg2.get_Molecule_by_name("PCT1")
PCT1.set_transition_width((0,1), width)
PCT2 = agg2.get_Molecule_by_name("PCT2")
PCT2.set_transition_width((0,1), width)
print("Aggregate has ", agg2.nmono, "single excited electronic states")

agg2.build(mult=2, vibgen_approx="TPA")

print("and ", agg2.Ntot, " (electro-vibrational) states in total")


# In[5]:


#
# Electronic aggregate is built with single exciton states only
#
agg_el.build(mult=1)
HHe = agg_el.get_Hamiltonian()


# In[6]:


#
#  Here we define system-bath interaction operators for relaxation in both the purely electronic
#  and the electro-vibrational Hamiltonian
#
#
e1_dim = agg_el.Nel
with qr.eigenbasis_of(HHe):

    K_12 = qr.qm.ProjectionOperator(1, 2, dim=e1_dim)
    K_21 = qr.qm.ProjectionOperator(2, 1, dim=e1_dim)
    K_23 = qr.qm.ProjectionOperator(2, 3, dim=e1_dim)
    K_32 = qr.qm.ProjectionOperator(3, 2, dim=e1_dim)
    K_24 = qr.qm.ProjectionOperator(2, 4, dim=e1_dim)
    #K1_toCT = qr.qm.ProjectionOperator(4, 1, dim=e1_dim)

    K_17 = qr.qm.ProjectionOperator(1, 5, dim=e1_dim)
    K_18 = qr.qm.ProjectionOperator(1, 6, dim=e1_dim)
    #print(K_12)

    # experimental kinetics
    ops = [K_12, K_23, K_24] #, K_17, K_18] #, K_21] #, K_21, K_32, K_13] #, K1_toCT]
    #rates = [1.0/100000000.0, 1.0/10000000.0, 1.0/100000] 
    kb_intK = qr.core.units.kB_intK
    expDEkbT = numpy.exp(-(HHe.data[2,2]-HHe.data[1,1])/(kb_intK*77.0))
    #print(expDEkbT,1.0/(expDEkbT*(1.0/27.0)))
    rates = [1.0/27.0, 1.0/157.0, 1.0/157.0] #, 1.0/300, 1.0/300] #, expDEkbT*(1.0/27.0)] #, 1.0/(27.0*20), 1.0/(157.0*20), 1.0/100.0] # 1.0/2100.0]

#
# Electronic only system-bath interaction operator
#
sbi_el = qr.qm.SystemBathInteraction(sys_operators=ops,
                                  rates=rates)
sbi_el.set_system(agg_el)

#
# System-bath interaction including vibrational states
#
sbi = qr.qm.SystemBathInteraction(sys_operators=ops,
                                  rates=rates)
sbi.set_system(agg)


# In[7]:


#
# Here we create Lindblad forms describing the relaxation
#

# Hamiltonian with vibrations (the purely electronic one was defined earlier)
HH = agg.get_Hamiltonian()

#
# Lindblad forms
#
LF_frac = qr.qm.ElectronicLindbladForm(HH, sbi, as_operators=True)
LF_el = qr.qm.ElectronicLindbladForm(HHe, sbi_el, as_operators=True)


# In[8]:


#
# We define time axis for propagation
# 
time_axis = qr.TimeAxis(0.0, 1000, 1.0)

#
# and propagator which can propagate purely electronic and electro-vibrational systems
#
prop = qr.qm.ReducedDensityMatrixPropagator(time_axis, HH, LF_frac)
prope = qr.qm.ReducedDensityMatrixPropagator(time_axis, HHe, LF_el)


# In[9]:


#
# Here we define the intial density matrix for propagation
#

# impulsive laser excitation, or single state population is the initial state
impulsive = True

if impulsive:

    #
    # Impulsive excitation using aggregate's transition dipole moments
    #
    rho0 = agg.get_DensityMatrix(condition_type="impulsive_excitation")
    rho0.normalize2()
    
    # we can trace out vibrations to get reduced density matrix which is electronic only
    # we check it against a purely electronic state below
    sig0 = agg.trace_over_vibrations(rho0)
    with qr.eigenbasis_of(HHe):
        # print populations
        print([numpy.real(sig0.data[i,i]) for i in range(1, e1_dim)])
    
    #
    # Impulsive excitation with purely electronic system
    #
    rho0e = agg_el.get_DensityMatrix(condition_type="impulsive_excitation")
    rho0e.normalize2()
    with qr.eigenbasis_of(HHe):
        # print populations
        print([numpy.real(rho0e.data[i,i]) for i in range(1, e1_dim)])
else:
 
    #
    # Single starting state
    #
    rho0 = agg.get_DensityMatrix(condition_type="impulsive_excitation")
    
    # we nullify the density matrix ...
    rho0.data[:,:] = 0.0
    with qr.eigenbasis_of(HH):
        # ... and set one state non-zero 
        rho0.data[6,6] = 1.0
    sig0 = agg.trace_over_vibrations(rho0)
 
    rho0e = agg_el.get_DensityMatrix(condition_type="impulsive_excitation")
    with qr.eigenbasis_of(HHe):    
        rho0e.data[:,:] = 0.0
        rho0e.data[2,2] = 1.0
    


# In[10]:


#
# THIS TAKES FEW SECONDS
#

#
# Propagation of the whole system including vibrations 
#
t1 = time.time()
rhot = prop.propagate(rho0)
t2 = time.time()
print("Propagation in:", t2-t1, "sec")

#
# Propagation of the electronic only system
#
t1 = time.time()
#with qr.eigenbasis_of(HHe):
rhoet = prope.propagate(rho0e)
t2 = time.time()
print("Propagation in:", t2-t1, "sec")


# In[11]:


#
# Time dependent density matrix with vibrations can be traced over the vibrations
#
sigt = agg.trace_over_vibrations(rhot)

# we check the shape of its data
print(rhot.data.shape)


# In[12]:


#
# Here we can compare electronic eigenstate populations at different times
#
Nt = 300
tm = time_axis.data[Nt]
print("Populations at time = ", tm)
with qr.eigenbasis_of(HHe):
    print("Traced: ")
    print([numpy.real(sigt.data[0,i,i]) for i in range(1, e1_dim)])
    print("Purely electronic: ")
    print([numpy.real(rhoet.data[0,i,i]) for i in range(1, e1_dim)])


# In[13]:


show_plots = True

#
# Plot of the dynamics
#
if show_plots:
    with qr.eigenbasis_of(HHe):
        #sigt.plot(coherences=False, axis=[0, 500, 0, 1.1], show=False)
        #rhoet.plot(coherences=False)
        plt.plot(rhot.TimeAxis.data, numpy.real(sigt.data[:,1,1]), "-r")
        plt.plot(rhot.TimeAxis.data, numpy.real(rhoet.data[:,1,1]), "-b")
        plt.show()
    
    with qr.eigenbasis_of(HH):
        plt.plot(rhot.TimeAxis.data, numpy.real(rhot.data[:,6,6]))
        plt.show()
    
    with qr.eigenbasis_of(HH):
        rhot.plot(coherences=False, show=True) 
    
    with qr.eigenbasis_of(HHe):
        rhoet.plot(coherences=False, show=True)    


# In[14]:


trD = agg_el.get_TransitionDipoleMoment()


# ## Liouville pathways

# In[15]:


#
# THIS TAKES FEW MINUTES (depending on N_steps)
#

#
#   Basis independent calculation of the evolution superoperator
#   The evolution superoperator is required for spectra calculations
#

N_steps = 10

#
# Time axis on which t2 is defined and FFT calculated
#
time_so = qr.TimeAxis(0.0, N_steps, 10.0)

eUt = qr.qm.EvolutionSuperOperator(time_so, HH, LF_frac)

#
# This cut the time step N times to make the numerics work !
#
eUt.set_dense_dt(10)

# This takes time (use eUt.calculate(show_progress=True) to see progress)
t1 = time.time()
eUt.calculate(show_progress=True)
t2 = time.time()
print("Finished in ", t2-t1, "sec")


# In[16]:


#
# Laboratory setup
#

from quantarhei import LabSetup
from quantarhei.utils.vectors import X, Y, Z

lab = LabSetup()
lab.set_polarizations(pulse_polarizations=[X,X,X], detection_polarization=X)

#
# Laser pulse weighting is on the way
#

#spectrum = qr.DFunction(time, values)
# lab.set_pulse_spectra((0,1,2,3), stype="power", spectrum )


# In[17]:


#
# THIS TAKES FEW TENS OF SECONDS
#
# In future version of Quantarhei, this call will not be needed (it will be done silently elsewhere, when needed)
t1 = time.time()
agg2.diagonalize()
t2 = time.time()
print("Diagonalized in ", t2-t1, "s")


# # Calculation of spectra at different $t_2$ times

# In[18]:


#
# we define a container for 2D spectra
#
cont = qr.TwoDSpectrumContainer(t2axis=time_so)
#
# spectra will be indexed by the times in the time axis `time_so`
#
cont.use_indexing_type(time_so)

#
# We define two-time axes, which will be FFTed and will define the omega_1 and omega_3 axes of the 2D spectrum
#
t1axis = qr.TimeAxis(0.0, 1000, 1.0)
t3axis = qr.TimeAxis(0.0, 1000, 1.0)

#
# This calculator calculated 2D spectra from the effective width defined above
#
msc = qr.MockTwoDSpectrumCalculator(t1axis, time_so, t3axis)
msc.bootstrap(rwa=qr.convert(12200.0,"1/cm","int"), 
              all_positive=False, shape="Gaussian")

#
# Now, let us calculate 2D spectra for all necessary time-points
#
tg1 = time.time()
for N_T2 in range(time_so.length):
    #
    #  Select t2 time
    #
    T2 = time_so.data[N_T2]
    print("\n***** Evolution at ", T2, " fs *****")
        
    #
    # Generation of Liouville pathways
    #
    print("Generating Liouville pathways")
    rho0 = agg2.get_DensityMatrix(condition_type="thermal", temperature=0.0)

    # types of pathways calculated
    reph = ("R2g", "R3g", "R1f*")
    noreph = ("R1g", "R4g", "R2f*")

    # add all types needed
    typs = reph #+ noreph
    
    #
    # Here we generate the pathways
    #
    t1 = time.time()
    pthways = agg2.liouville_pathways_3T(ptype=typs,
                                          lab=lab,
                                          eUt=eUt, t2=T2, etol=1.0e-5, verbose=0) #eUt2=qr.qm.SOpUnity(dim=HH.dim))
    t2 = time.time()
    print("Generation time: ", t2-t1, "sec")
    print("Number of pathways: ", len(pthways))

    #
    # Let's process the pathways
    #
    t1 = time.time()
    print("Selecting relevant pathways")
    pw = []
    
    # we take a range of t2 frequencies
    om_low = qr.convert(570, "1/cm", "int")
    om_up = qr.convert(575, "1/cm", "int")
    pw_p570 = spec.select_omega2((om_low, om_up), pthways)
    pw_m570 = spec.select_omega2((-om_up, -om_low), pthways)

    # take only pathways with amplitude larger than a value
    pw_p = spec.select_amplitude_GT(0.000001, pw_p570)
    pw_m = spec.select_amplitude_GT(0.000001, pw_m570)
    pw_S = []
    pw_S += pw_p
    pw_S += pw_m
    pw = spec.order_by_amplitude(pw_S)

    #
    # Select pathways from different spectral regions
    #
    
    # on P-
    pw_Pm = spec.select_frequency_window([qr.convert(11000, "1/cm", "int"), qr.convert(11500, "1/cm", "int"),
                                   qr.convert(11000, "1/cm", "int"), qr.convert(13800, "1/cm", "int")], pw)

    # on P+
    pw_Pp = spec.select_frequency_window([qr.convert(11500, "1/cm", "int"), qr.convert(12200, "1/cm", "int"),
                                   qr.convert(11000, "1/cm", "int"), qr.convert(13800, "1/cm", "int")], pw)

    # on B
    pw_B = spec.select_frequency_window([qr.convert(12200, "1/cm", "int"), qr.convert(13000, "1/cm", "int"),
                                   qr.convert(11000, "1/cm", "int"), qr.convert(13800, "1/cm", "int")], pw)

    print("Number of pathways on P- ", len(pw_Pm))
    print("Number of pathways on P+ ", len(pw_Pp))
    print("Number of pathways on B  ", len(pw_B))

    #
    # Which pathways to present
    #
    pw = []
    pw += pw_Pm
    pw += pw_Pp
    pw += pw_B

    t2 = time.time()
    print("Pathways selected in", t2-t1, "sec")
    
    t1 = time.time()
    print("Calculating frequency map")

    #
    # save pathways
    #
    pw.save("pathways_"+str(T2)+".h5py")

    #
    # Calculate 2D spectrum using the predefined calculator
    #
    msc.set_pathways(pw)
    twod = msc.calculate_next()
    t2 = time.time()
    print("..done in", t2-t1,"sec")
    
    #
    # set current t2 time to the spectrum as a tag
    #
    twod.set_t2(T2)
    
    #
    # And put the spectrum to the container
    #
    cont.set_spectrum(twod)

tg2 = time.time()
print("In total it took: ", tg2-tg1, "sec.")


# In[19]:


#
# Here you can inspect and/or recalculate some spectra
#
# If repeat is True, you can construct a new spectrum from existing pathways,
# i.e. those regarding the last calculated spectrum
#
#

repeat = False

if repeat:

    pw_p = spec.select_amplitude_GT(0.000001, pw_p570)
    pw_m = spec.select_amplitude_GT(0.000001, pw_m570)
    pw_S = []
    pw_S += pw_p
    pw_S += pw_m
    pw = spec.order_by_amplitude(pw_S)

    # on P-
    pw_Pm = spec.select_frequency_window([qr.convert(11000, "1/cm", "int"), qr.convert(11500, "1/cm", "int"),
                                   qr.convert(11000, "1/cm", "int"), qr.convert(13800, "1/cm", "int")], pw)

    # on P+
    pw_Pp = spec.select_frequency_window([qr.convert(11500, "1/cm", "int"), qr.convert(12200, "1/cm", "int"),
                                   qr.convert(12700, "1/cm", "int"), qr.convert(13800, "1/cm", "int")], pw)

    # on B
    pw_B = spec.select_frequency_window([qr.convert(12200, "1/cm", "int"), qr.convert(13000, "1/cm", "int"),
                                   qr.convert(11000, "1/cm", "int"), qr.convert(13800, "1/cm", "int")], pw)

    print("Number of pathways on P- ", len(pw_Pm))
    print("Number of pathways on P+ ", len(pw_Pp))
    print("Number of pathways on B  ", len(pw_B))    
    
    pw = []
    #pw += pw_Pm
    pw += pw_Pp
    #pw += pw_B

    t2 = time.time()
    print("Pathways selected in", t2-t1, "sec")
    
    t1 = time.time()
    print("Calculating frequency map")

    msc.bootstrap(rwa=qr.convert(12200.0,"1/cm","int"), 
              all_positive=False, shape="Gaussian")
    msc.set_pathways(pw)
    twod = msc.calculate_next()
    t2 = time.time()
    print("..done in", t2-t1,"sec")
    
else:
    
    twod = cont.get_spectrum(10.0)



#
# When this cell is left, we either retrieve one spectrum from container, or calculate a new one from available pathways
#


# In[20]:


#
# Plotting an example spectrum
#
with qr.energy_units("1/cm"):
    twod.plot(window=[11000,14000,11000,14000], Npos_contours=10, 
              stype="total", spart="real")


# In[21]:


#
# If the plotting window is reasonable, we should cut the unnecessary data by trimming the 2D spectrum
#
with qr.energy_units("1/cm"): 
    cont.trimall_to(window=[11000,14000,11000,14000])


# In[22]:


#
# twod spectrum object is under revision - use d__data instead of data if you want to reach the raw data
#
sp = cont.get_spectrum(10.0)
sp.set_data_flag("total")
print(sp.d__data[33,30])


# In[23]:


#
# Window function for subsequenty FFT
#
import quantarhei.functions as func
window = func.Tukey(time_so, r=0.3, sym=False)

#
# FFT with the window function
#
# Specify REPH, NONR or `total` to get different types of spectra
#
fcont = cont.fft(window=window, dtype="REPH")


# In[24]:


#
#  Here we save all results
#

# save container with spectra
cont.save("all_spectra.h5py")
fcont.save("all_spectra_fft.h5py")
agg2.save("aggregate.h5py")
eUt.save("evol_operator.h5py")





#
# Have a look which frequencies we actually have
#
Ndat = len(fcont.axis.data)
print("Number of frequency points:", Ndat)
print("\nIn 1/cm they are:")
with qr.energy_units("1/cm"):
    for k_i in range(Ndat):
        print(k_i, fcont.axis.data[k_i])
    


# In[25]:


#
# Current storage of the specta in the container works on string representation of the frequency. 
# As a consequence the container does not recognize units of frequency. We print the frequency 
# in 1/cm, but we have to specify it in internal units
#

Npoint = 3
om = fcont.axis.data[Npoint]
sp = fcont.get_spectrum(om)
units = "1/cm"
with qr.energy_units(units):
    print("Spectrum at frequency:", fcont.axis.data[Npoint], units)
    sp.plot(window=[11000,14000,11000,14000], Npos_contours=20, 
              stype="total", spart="abs")


# # Check individual pathways here

# In[26]:


with qr.energy_units("1/cm"):
    for pathw in pw:
        print(pathw)


# In[27]:


with qr.energy_units("1/cm"):
    pathw = pw[0]
    print(pathw)


# In[28]:


agg2.report_on_expansion(10)


# In[29]:


agg2.report_on_expansion(40)


# In[30]:


N1 = 40
N2 = 37
with qr.energy_units("1/cm"):
    print(agg2.get_state_energy(N1))
    print(agg2.get_state_energy(N2))
    print(agg2.get_transition_dipole(N1))
    print(agg2.get_transition_dipole(N2))


# In[31]:


agg2.report_on_expansion(123)


# In[32]:


agg2.report_on_expansion(36)


# In[33]:


agg2.report_on_expansion(8)


# In[34]:


agg2.report_on_expansion(6)


# In[35]:


agg2.report_on_expansion(27)


# In[36]:


agg2.report_on_expansion(28)


# In[37]:


agg2.report_on_expansion(3)


# In[38]:


agg2.report_on_expansion(50)


# In[39]:


agg2.report_on_expansion(18)


# In[40]:


HH2 = agg2.get_Hamiltonian()


# In[41]:


with qr.energy_units("1/cm"):
    print([HH2.data[i,i] for i in range(4)])


# # Time evolution of the components

# In[42]:



plt.plot(time_so.data, numpy.real(eUt.data[:,40,37,10,7]), "-g")

