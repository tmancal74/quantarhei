# coding: utf-8

#
# Reading configuration file 
#
#
with open("fft2D.conf","r") as f:
    txt = f.read()
    exec(txt)
    
# try the configuration
import traceback
import datetime
try:
    print("\nRC Simulation", datetime.datetime.now(),"\n")
    print("Simulation parameters: ")
    print("restart  = ", restart)
    print("eUt_mode = ", eUt_mode)
    
    
    print("propagation_dt = ", propagation_dt)
    print("Using pure dephasing = ", pure_deph)
    print("")
    
    print("2D FFT parameters:")
    print("FFT offset: ", t_offset)
    print("Plot window (in cm-1) = ", plot_window)
    print("")    
except:
    traceback.print_exc()
    print("Configuration file is incomplete")
    quit()


# In[1]:

import os

import quantarhei as qr
import quantarhei.spectroscopy as spec
import numpy
import time

import matplotlib.pyplot as plt
plt.switch_backend('agg')

print("Quantarhei version: ", qr.Manager().version, "\n")

# input and output paths of the script
pre_in = os.path.join("..", "model", "out")
pre_out = "out"

# check if pre_out exists and is a directory
if not os.path.isdir(pre_out):
    try:
        os.makedirs(pre_out, exist_ok=True)
    except:
        raise Exception("Output directory name '"
                        +pre_out+"' does not represent a valid directory")

# In[2]:


#
# This is needed in version 0.0.36 for propagation with Lindblad form
#
qr.Manager().gen_conf.legacy_relaxation = True



# In[4]:

#
# Load aggregates constructed in the previous phase of the calculation
#

agg = qr.load_parcel(os.path.join(pre_in, 
                                  "fraction_45_2_vibrations_CT_unbuilt.qrp"))
agg2 = qr.load_parcel(os.path.join(pre_in,
                                   "fraction_45_2_vibrations_CT_unbuilt.qrp"))
agg_el = qr.load_parcel(os.path.join(pre_in,"fraction_eff_40_4_CT_unbuilt.qrp"))   


# In[5]:


#
# Aggregate with vibrational states, Hamiltonian generated up to single
# exciton states in a Two-particle approximation
# This object is used for calculation of dynamics
#
agg.build(mult=1, vibgen_approx="TPA")

#
# agg2 object will be used for calculation of 2D spectra
# so it needs to know how to represent the spectra (what width it should have).
# It also has to be built with
# two-exciton states (in two-particle approximation)
#
width = qr.convert(wincm_P, "1/cm", "int")
PM = agg2.get_Molecule_by_name("PM")
PM.set_transition_width((0,1), width)
PL = agg2.get_Molecule_by_name("PL")
PL.set_transition_width((0,1), width)

width = qr.convert(wincm_B, "1/cm", "int")
BM = agg2.get_Molecule_by_name("BM")
BM.set_transition_width((0,1), width)
BL = agg2.get_Molecule_by_name("BL")
BL.set_transition_width((0,1), width)

width = qr.convert(wincm_CT, "1/cm", "int")
PCT1 = agg2.get_Molecule_by_name("PCT1")
PCT1.set_transition_width((0,1), width)
PCT2 = agg2.get_Molecule_by_name("PCT2")
PCT2.set_transition_width((0,1), width)
print("Aggregate has ", agg2.nmono, "single excited electronic states")

agg2.build(mult=2, vibgen_approx="TPA")

print("and ", agg2.Ntot, " (electro-vibrational) states in total")


# In[6]:


#
# Electronic aggregate is built with single exciton states only
#
width = qr.convert(wincm_P, "1/cm", "int")
PM = agg_el.get_Molecule_by_name("PM")
PM.set_transition_width((0,1), width)
PL = agg_el.get_Molecule_by_name("PL")
PL.set_transition_width((0,1), width)

width = qr.convert(wincm_B, "1/cm", "int")
BM = agg_el.get_Molecule_by_name("BM")
BM.set_transition_width((0,1), width)
BL = agg_el.get_Molecule_by_name("BL")
BL.set_transition_width((0,1), width)

width = qr.convert(wincm_CT, "1/cm", "int")
PCT1 = agg_el.get_Molecule_by_name("PCT1")
PCT1.set_transition_width((0,1), width)
PCT2 = agg_el.get_Molecule_by_name("PCT2")
PCT2.set_transition_width((0,1), width)
print("Aggregate has ", agg_el.nmono, "single excited electronic states")
agg_el.build(mult=1)
HHe = agg_el.get_Hamiltonian()


# In[7]:


#
#  Here we define system-bath interaction operators for relaxation in both
#  the purely electronic and the electro-vibrational Hamiltonian
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

    # thermal factor for the rates
    kb_intK = qr.core.units.kB_intK
    expDEkbT = numpy.exp(-(HHe.data[2,2]-HHe.data[1,1])/(kb_intK*77.0))

    # experimental kinetics
    ops = [K_12, K_23, K_24] 
        
    # rates in the order defined above
    rates = [1.0/27.0, 1.0/157.0, 1.0/157.0] 

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


# In[8]:


#
# Here we create Lindblad forms describing the relaxation
#

# Hamiltonian with vibrations (the purely electronic one was defined earlier)
HH = agg.get_Hamiltonian()

#
# Lindblad forms
#

# with vibrations
LF_frac = qr.qm.ElectronicLindbladForm(HH, sbi, as_operators=True)

# electronic only
LF_el = qr.qm.ElectronicLindbladForm(HHe, sbi_el, as_operators=True)


# In[9]:
#
# Laboratory setup
#

from quantarhei import LabSetup
from quantarhei.utils.vectors import X, Y, Z

lab = LabSetup()
lab.set_polarizations(pulse_polarizations=[X,X,X], detection_polarization=X)

#
# We define time axis for propagation
# 
time_axis = qr.TimeAxis(0.0, propagation_N_steps, propagation_dt)

agg_el.diagonalize()

#
#  Absorption Spectrum by pathway method
#

mac = qr.MockAbsSpectrumCalculator(time_axis, system=agg_el)

rho0 = agg_el.get_DensityMatrix(condition_type="thermal", temperature=0.0)
ham = agg_el.get_Hamiltonian()

pthways = agg_el.liouville_pathways_1(lab=lab, ham=ham, etol=1.0e-5,
                                       verbose=0) 

mac.bootstrap(rwa=qr.convert(12200.0,"1/cm","int"), 
              shape="Gaussian")

mac.set_pathways(pthways)

abs1 = mac.calculate()
abs1.normalize2(norm=0.53)

#
# Absorption by standard method
#
try:
    abscalc = qr.load_parcel("../model/in/calcAbs.qrp")
except:
    pass
#
# Plotting absorption spectra
#
try:
    absexp = qr.load_parcel("../model/in/bas_77K.qrp")

    absexp.normalize()
    absexp.subtract(0.086)
    absexp.normalize2(norm=0.53)
except:
    pass

plt.figure(1)
with qr.energy_units("1/cm"):
    abs1.plot(axis=[10500,15000, 0.0, 0.7], show=False)
    try:
        absexp.plot()
    except:
        pass
    try:
        abscalc.plot()
    except:
        pass
    
plt.savefig(os.path.join(pre_out,"abs1.png"))
plt.close()



#
# and propagator which can propagate purely electronic and electro-vibrational systems
#
prop = qr.qm.ReducedDensityMatrixPropagator(time_axis, HH, RTensor=LF_frac)
prope = qr.qm.ReducedDensityMatrixPropagator(time_axis, HHe, RTensor=LF_el)


# In[10]:


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
    
    # we can trace out vibrations to get reduced density matrix which
    # is electronic only
    # we check it against a purely electronic state below
    print("Comparison of pure electronic and vibrational but traced over DMs")
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
    



# In[11]:


#
# THIS TAKES FEW SECONDS
#

print("Propagating the system (electronic only and vibrational)")

#
# Propagation of the whole system including vibrations 
#
print("Including vibrations ...")
t1 = time.time()
rhot = prop.propagate(rho0)
t2 = time.time()
print("Propagation in:", t2-t1, "sec")

#
# Propagation of the electronic only system
#
print("Electronic only ...")
t1 = time.time()
#with qr.eigenbasis_of(HHe):
rhoet = prope.propagate(rho0e)
t2 = time.time()
print("Propagation in:", t2-t1, "sec")


# In[12]:


#
# Time dependent density matrix with vibrations can be traced over the vibrations
#
sigt = agg.trace_over_vibrations(rhot)

# we check the shape of its data
#print(rhot.data.shape)


# In[13]:


#
# Here we can compare electronic eigenstate populations at different times
#
print("Checking populations at one point:")
Nt = int(propagation_N_steps/2)
tm = time_axis.data[Nt]
print("Populations at time = ", tm)
with qr.eigenbasis_of(HHe):
    print("Traced: ")
    print([numpy.real(sigt.data[0,i,i]) for i in range(1, e1_dim)])
    print("Purely electronic: ")
    print([numpy.real(rhoet.data[0,i,i]) for i in range(1, e1_dim)])


# In[14]:


show_plots = True

#
# Plot of the dynamics
#
if show_plots:
    with qr.eigenbasis_of(HHe):
        #sigt.plot(coherences=False, axis=[0, 500, 0, 1.1], show=False)
        #rhoet.plot(coherences=False)
        plt.figure(1)
        plt.plot(rhot.TimeAxis.data, numpy.real(sigt.data[:,1,1]), "-r")
        plt.plot(rhot.TimeAxis.data, numpy.real(rhoet.data[:,1,1]), "-b")
        #plt.show()
        plt.savefig(os.path.join(pre_out,"fig01_kinetics_electronic_1_1.png"))
        plt.close()
    
    with qr.eigenbasis_of(HH):
        plt.figure(1)
        plt.plot(rhot.TimeAxis.data, numpy.real(rhot.data[:,6,6]))
        #plt.show()
        plt.savefig(os.path.join(pre_out,"fig02_kinetics_all_6_6.png"))
        plt.close()
        
    with qr.eigenbasis_of(HH):
        plt.figure(1)
        rhot.plot(coherences=False, show=False)
        #plt.show()
        plt.savefig(os.path.join(pre_out,"fig03_kinetics_all_all.png"))
        plt.close()
        
    with qr.eigenbasis_of(HHe):
        plt.figure(1)
        rhoet.plot(coherences=False, show=False)
        #plt.show()
        plt.savefig(os.path.join(pre_out,"fig04_kinetics_electronic_all.png"))
        plt.close()


# In[15]:


trD = agg_el.get_TransitionDipoleMoment()


# ## Liouville pathways

# In[16]:

if restart:
    
    print("\nReloading last state ")
    try:
        state = qr.load_parcel("A_saved_state.qrp")
        N_T2 = state[0]  # here we take the saved value
        agg2 = state[3]
        print("... reload succesful")
    except:
        N_T2 = 0
        restart = False
        print("... reload failed; starting from the beginning")

else:
    N_T2 = 0

#
# THIS TAKES FEW TENS OF SECONDS
#
# In future version of Quantarhei, this call will not be needed (it will be done silently elsewhere, when needed)
print("\nDiagonalization of the aggregate representation:")
t1 = time.time()
agg2.diagonalize()
t2 = time.time()
print("Diagonalized in ", t2-t1, "s")


#pure_deph = True
#t_offset = 80.0

if pure_deph:
    print("Calculating electronic pure dephasing")
    t1 = time.time()
    p_deph = qr.qm.ElectronicPureDephasing(agg2, dtype="Gaussian")
    t2 = time.time()
    print("Pure dephasing calculated in ", t2-t1, "s")
else:
    p_deph = None

#
# THIS TAKES FEW MINUTES (depending on t2_N_steps)
#

#
#   Basis independent calculation of the evolution superoperator
#   The evolution superoperator is required for spectra calculations
#

print("Setting up evolution superoperator at", t2_N_steps, "points:")
#
# Time axis on which t2 is defined and FFT calculated
#
time_so = qr.TimeAxis(0.0, t2_N_steps, t2_time_step)

eUt = qr.qm.EvolutionSuperOperator(time_so, HH, relt=LF_frac, 
                                   pdeph=p_deph, mode=eUt_mode)

#
# This cuts the time step N times to make the numerics work !
#
eUt.set_dense_dt(t2_sub)

#eUt_mode = "jit" # THIS IS SET IN A CONFIGURATION FILE

if eUt_mode == "all":
    # This takes time (use eUt.calculate(show_progress=True) to see progress)
    print("Calculating the whole dynamics in advance ...")
    t1 = time.time()
    eUt.calculate(show_progress=True)
    t2 = time.time()
            
#    indx = (eUt.data[2,:,:,:,:]>1.0).nonzero()
#    no_elem = len(indx[0])
#    for ii in range(no_elem):
#        print("*** - ", (indx[0][ii], indx[1][ii], indx[2][ii], indx[3][ii]))
#    eUt.plot_element((indx[0][0], indx[1][0], indx[2][0], indx[3][0]))
    eUt.plot_element((0, 3, 0, 9))
    eUt.plot_element((10, 10, 13, 13))
    eUt.plot_element((10, 13, 10, 13))
    plt.savefig(os.path.join(pre_out,"element.png"))
    
    print("Finished in ", t2-t1, "sec")
else:
    print("Dynamics will be calculated on fly")

# In[17]:

print("Calculating 2D spectra:")


#
# Laser pulse weighting is on the way
#

#spectrum = qr.DFunction(time, values)
# lab.set_pulse_spectra((0,1,2,3), stype="power", spectrum )


# In[18]:



# # Calculation of spectra at different $t_2$ times

# In[19]:



#
# we define a container for 2D spectra
#
cont = qr.TwoDSpectrumContainer(t2axis=time_so)
#
# spectra will be indexed by the times in the time axis `time_so`
#
cont.use_indexing_type(time_so)

#
# We define two-time axes, which will be FFTed and will define the omega_1 and
# omega_3 axes of the 2D spectrum
#
t1axis = qr.TimeAxis(0.0, t1_N_steps, t1_time_step)
t3axis = qr.TimeAxis(0.0, t3_N_steps, t3_time_step)

#
# This calculator calculated 2D spectra from the effective width defined above
#
msc = qr.MockTwoDSpectrumCalculator(t1axis, time_so, t3axis)
msc.bootstrap(rwa=qr.convert(12200.0,"1/cm","int"), 
              all_positive=False, shape="Gaussian")

ham = eUt.get_Hamiltonian()
    
    
#
# Now, let us calculate 2D spectra for all necessary time-points
#
print("2D spectra ...")
tg1 = time.time()

while (N_T2 < time_so.length):
#for N_T2 in range(time_so.length):
    
    #
    #  Select t2 time
    #
    T2 = time_so.data[N_T2]
    print("\n***** Evolution at ", T2, " fs *****")
        
          
    if not restart:

        #
        # Generation of Liouville pathways
        #
    
        rho0 = agg2.get_DensityMatrix(condition_type="thermal", temperature=0.0)
    
        # types of pathways calculated
        reph = ("R2g", "R3g", "R1f*")
        noreph = ("R1g", "R4g", "R2f*")
    
        # add all types needed
        typs = reph #+ noreph
            
        #
        # Here we generate the pathways
        #
        print("Generating Liouville pathways ...")    
        t1 = time.time()
        if eUt_mode == "all":
            eUt2 = eUt.at(T2)
        else:
            eUt2 = eUt.at()
        
        pthways = agg2.liouville_pathways_3T(ptype=typs,
                                              lab=lab,
                                              eUt=eUt2, ham=ham, t2=T2, 
                                              etol=1.0e-6,
                                              verbose=0) 
                                              #eUt2=qr.qm.SOpUnity(dim=HH.dim))
        t2 = time.time()
        print(" ... done")
        print("Generation time: ", t2-t1, "sec")
        print("Number of pathways: ", len(pthways))
    
        #
        # Let's process the pathways
        #
        t1 = time.time()
        print("Selecting relevant pathways")
        pw = []
        
        # we take a range of t2 frequencies
        om_low = qr.convert(550, "1/cm", "int")
        om_up = qr.convert(650, "1/cm", "int")
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
        pw_Pm = spec.select_frequency_window(
                [qr.convert(11000, "1/cm", "int"), 
                 qr.convert(11500, "1/cm", "int"),
                 qr.convert(11000, "1/cm", "int"), 
                 qr.convert(13800, "1/cm", "int")], pw)
    
        # on P+
        pw_Pp = spec.select_frequency_window(
                [qr.convert(11500,"1/cm", "int"),
                 qr.convert(12200, "1/cm", "int"),
                 qr.convert(11000, "1/cm", "int"),
                 qr.convert(13800, "1/cm", "int")], pw)
    
        # on B
        pw_B = spec.select_frequency_window(
                [qr.convert(12200, "1/cm", "int"),
                 qr.convert(13000, "1/cm", "int"),
                 qr.convert(11000, "1/cm", "int"),
                 qr.convert(13800, "1/cm", "int")], pw)

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
        # Calculate 2D spectrum using the predefined calculator
        #
        msc.set_pathways(pw)
        twod = msc.calculate_next()
        t2 = time.time()
        print("... done in", t2-t1,"sec")
        
        #
        # set current t2 time to the spectrum as a tag
        #
        twod.set_t2(T2)
        
        #
        # And put the spectrum to the container
        #
    
        cont.set_spectrum(twod)
    
        print("Saving all pathways:")
        prcl = qr.save_parcel(pw,
                              os.path.join(pre_out,
                                           "pathways"+"_"+str(T2)+".qrp"),
                              comment="pathways at t2 = "+str(T2))
        print("... done")
    
    
        #
        # Save stuff for a restart
        #
        if eUt_mode == "jit":
            state = [N_T2, cont, eUt, agg2] #, msc]
            prcl = qr.Parcel()
            prcl.set_content(state)
            prcl.save("A_saved_state.qrp")
            
    else:
        
        cont = state[1]
        eUt = state[2]
        #msc = state[4]
        if N_T2+1 < time_so.length:
            print("Setting next step as ", N_T2+1, " i.e. ",
                  time_so.data[N_T2+1], "fs")
            msc.set_next(N_T2+1) 
        
        cont.t2axis = time_so
        cont.use_indexing_type(time_so)
        restart = False
        

    #
    # propagation of the evolution superoperator
    #
    if eUt_mode == "jit":
        print("Propagating dynamics ...")
        t1 = time.time()
    
        if eUt.has_PureDephasing():
            # pure dephasing has to be used in site basis
            print("   calculating in eigenstate basis")
            with qr.eigenbasis_of(ham):
                eUt.calculate_next()
        else:
            # here we avoid unnecessary basis transformations
            eUt.calculate_next()
        
        t2 = time.time()
        print("... propagated in ", t2-t1, "sec")
        
    N_T2 += 1

tg2 = time.time()
print("In total 2D spectra calculation took: ", tg2-tg1, "sec.")




# In[20]:


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
    pw_Pm = spec.select_frequency_window([qr.convert(11000, "1/cm", "int"),
                                          qr.convert(11500, "1/cm", "int"),
                                          qr.convert(11000, "1/cm", "int"),
                                          qr.convert(13800, "1/cm", "int")], pw)

    # on P+
    pw_Pp = spec.select_frequency_window([qr.convert(11500, "1/cm", "int"),
                                          qr.convert(12200, "1/cm", "int"),
                                          qr.convert(12700, "1/cm", "int"),
                                          qr.convert(13800, "1/cm", "int")], pw)

    # on B
    pw_B = spec.select_frequency_window([qr.convert(12200, "1/cm", "int"),
                                         qr.convert(13000, "1/cm", "int"),
                                         qr.convert(11000, "1/cm", "int"),
                                         qr.convert(13800, "1/cm", "int")], pw)

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
    print("... done in", t2-t1,"sec")
    
else:
    
    N_T2_pul = int(N_T2/2)
    print("Saving 2D spectrum at",
          time_so.data[N_T2_pul], "fs (as an example)")
    twod = cont.get_spectrum(time_so.data[N_T2_pul])

#
# When this cell is left, we either retrieve one spectrum from container, or calculate a new one from available pathways
#


# In[21]:


#
# Plotting an example spectrum
#
    
# plot_window = [11000,13500,11000,13500]
    
with qr.energy_units("1/cm"):
    #plt.figure(1)
    twod.plot(window=plot_window, Npos_contours=10,              
              stype="total", spart="real")
    plt.savefig(os.path.join(pre_out,"twod_example_at_T2.png"))
    #plt.close()


# In[22]:


#
# If the plotting window is reasonable, we should cut the unnecessary data by trimming the 2D spectrum
#
with qr.energy_units("1/cm"): 
    cont.trimall_to(window=plot_window)


# In[24]:


#
# Window function for subsequenty FFT
#
import quantarhei.functions as func
window = func.Tukey(time_so, r=tukey_r, sym=False)

#
# FFT with the window function
#
# Specify REPH, NONR or `total` to get different types of spectra
#
print("Calculating FFT of the 2D maps")
fcont = cont.fft(window=window, dtype="REPH", offset=t_offset)


# In[25]:


#
# Have a look which frequencies we actually have
#
Ndat = len(fcont.axis.data)
print("Number of frequency points:", Ndat)
print("\nIn 1/cm they are:")
with qr.energy_units("1/cm"):
    for k_i in range(Ndat):
        print(k_i, fcont.axis.data[k_i])
    


# In[26]:


#
# Current storage of the specta in the container works on string representation of the frequency. 
# As a consequence the container does not recognize units of frequency. We print the frequency 
# in 1/cm, but we have to specify it in internal units
#

# The point with the frequency nearest to the desired one should be chosen
#show_Npoint = 7 #int(3*Ndat/4)

om = fcont.axis.data[show_Npoint]
#print("omega: ", om)
sp = fcont.get_spectrum(om)

units = "1/cm"
with qr.energy_units(units):
    print("Spectrum at frequency:", fcont.axis.data[show_Npoint], units)
    sp.plot(window=plot_window, Npos_contours=30, 
              stype="total", spart="abs")
    sp.savefig(os.path.join(pre_out, "twod_fft_map.png"))

# In[46]:

print("Saving spectral containers - both time and fft")
qr.save_parcel(fcont,os.path.join(pre_out,"spectra_container_fft.qrp"))
qr.save_parcel(cont,os.path.join(pre_out,"spectra_container_time.qrp"))


# # Check individual pathways here

# In[28]:


skip = True

if not skip:
    with qr.energy_units("1/cm"):
        for pathw in pw:
            print(pathw)


# In[29]:

    print("The most intensive pathway (as an example)")
    with qr.energy_units("1/cm"):
        pathw = pw[0]
        print(pathw)


#print("Saving evolution superoperator")
#peUt = qr.Parcel()
#peUt.set_content(eUt)
#peUt.save(os.path.join(pre_out,"evolution_op_eUt.qrp"))

