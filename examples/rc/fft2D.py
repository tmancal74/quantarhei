# coding: utf-8
###############################################################################
#
#
#
#
#
#
#
###############################################################################

#
# Imports
#

# standard libraries
import os
import datetime
import time

# pyplot
import matplotlib.pyplot as plt
plt.switch_backend('agg')

try:
    from terminaltables import AsciiTable
except:
    raise Exception("Get terminaltables package "
                            +"for this functionality")

# numerics
import numpy

# quantarhei
import quantarhei as qr
import quantarhei.spectroscopy as spec

#
###############################################################################
#
# SCRIPT CONFIGURATION
#
###############################################################################
#


#
# Reading of configuration file 
#
with open("fft2D.conf","r") as f:
    txt = f.read()
    exec(txt)
 
def tprint(var, messg=None, default=None):
    """Test the existence of the variable with the name `var`
    
    
    If the default is specified, the non-existence of the variable
    is not a problem, and the variable is set to the default
    
    """    
    import traceback
    
    try:
        if messg is not None:
            print("#", messg)
        val = eval(var)
        if isinstance(val, str):
            val = '"'+val+'"'
        print(var, "=", val)

    except:
        if default is None:
            traceback.print_exc()
            print("Configuration file is incomplete")
            quit()
        else:
            globals()[var] = default
            val = eval(var)
            if isinstance(val, str):
                val = '"'+val+'"'
            print(var,"=", val, "# (default)")

    
#
# Test the compulsory part of the configuration
#

print("")     
print("######################################################################")
print("#                                                                    #")
print("#             FT VIBRONIC 2D MAPS of REACTION CENTER                 #")
print("#                                                                    #")
print("#             Quantarhei simulation script                           #")
print("#                                                                    #")
print("#    This script is a part of the supporting information             #")
print("#    of the following publication:                                   #")
print("#                                                                    #")      
print("#    Policht, V.; ..., (2019)                                        #")
print("#                                                                    #")
print("#                                                                    #")      
print("######################################################################")
print("")
print("Simulation started at: "+str(datetime.datetime.now()))
print("Quantarhei version: ", qr.Manager().version)
print("")
print("***                      Simulation parameters                     ***")
print("\n# This block can be used as an input file for another calculation")
print("\n# Input and output directories")
tprint("out_dir", default="out")
tprint("model", default=os.path.join("..","model"))
print("\n# Waiting time propagation parameters:")
tprint("eUt_mode", default="jit")
tprint("restart", default=True)
tprint("stop_after_propagation", default=False)
tprint("pure_deph", messg="Using pure dephasing?", default=True)
tprint("t2_N_steps")
tprint("t2_time_step")
tprint("t2_sub", default=10)
tprint("transfer_times")
tprint("uphill_thermal", default = True)
tprint("ignore_slower_than", messg="Ignore uphill rates slower than this (fs)")

print("\n# Time domain 2D calculation parameters:")
tprint("t1_N_steps", messg="t1 coherence time")
tprint("t1_time_step")
tprint("t3_N_steps", messg="t3 coherence time")
tprint("t3_time_step")

print("\n# Effective lineshape parameters:")
tprint("transition_widths")

print("\n# Demo propagation parameters:")
tprint("propagation_dt")
tprint("propagation_N_steps")

print("\n# 2D FFT parameters:")
tprint("plot_window")
tprint("t_offset", default=0.0, messg="FFT offset")
tprint("tukey_r")
tprint("show_omega", default=0.0)
print("\n# Liouville pathway selection:")
tprint("frequency_interval")
if show_omega < frequency_interval[0] or show_omega > frequency_interval[1]:
    print("  WARNING: pathways will not have right frequencies (show_omega)"+
          "\n  is out of selection window")
    raise Exception("Stopping because of a warning")

print("\n# Auxiliary plots:")
tprint("evol_super_op_elems2plot", default=[],
       messg="Evolution Superoperator Elements to be plotted")
tprint("show_kinetics", default=False)
tprint("comparison_elems", default=(1,1))
tprint("detail_elems", default=(6,6))

print("\n# End of the input file")
print("")
print("***                      Simulation output                         ***")
print("")

pre_out = out_dir
pre_in = os.path.join(model, "out")


# check if pre_out exists and is a directory
if not os.path.isdir(pre_out):
    try:
        os.makedirs(pre_out, exist_ok=True)
    except:
        raise Exception("Output directory name '"
                        +pre_out+"' does not represent a valid directory")


#
# This is needed in version 0.0.36 and later for propagation with Lindblad form
#
qr.Manager().gen_conf.legacy_relaxation = True


#
# Load aggregates constructed in the previous phase of the calculation
#

agg = qr.load_parcel(os.path.join(pre_in, 
                                  "fraction_45_2_vibrations_CT_unbuilt.qrp"))
agg2 = qr.load_parcel(os.path.join(pre_in,
                                   "fraction_45_2_vibrations_CT_unbuilt.qrp"))
agg_el = qr.load_parcel(os.path.join(pre_in,
                                     "fraction_eff_40_4_CT_unbuilt.qrp"))   

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

for mol_name in transition_widths:
    
    width = qr.convert(transition_widths[mol_name], "1/cm", "int")
    mol = agg2.get_Molecule_by_name(mol_name)
    mol.set_transition_width((0,1), width)

print("Aggregate has ", agg2.nmono, "single excited electronic states")

agg2.build(mult=2, vibgen_approx="TPA")

print("and ", agg2.Ntot, " (electro-vibrational) states in total")
print("Number of single exciton states :", agg2.Nb[0]+agg2.Nb[1])


qr.save_parcel(agg2, os.path.join(pre_out,"agg2_built.qrp"))

#
# Electronic aggregate is built with single exciton states only
#

for mol_name in transition_widths:
    
    width = qr.convert(transition_widths[mol_name], "1/cm", "int")
    mol = agg_el.get_Molecule_by_name(mol_name)
    mol.set_transition_width((0,1), width)

print("Aggregate has ", agg_el.nmono, "single excited electronic states")
agg_el.build(mult=1)
HHe = agg_el.get_Hamiltonian()

#
#  Here we define system-bath interaction operators for relaxation in both
#  the purely electronic and the electro-vibrational Hamiltonian
#
#

kb_intK = qr.core.units.kB_intK
e1_dim = agg_el.Nel
with qr.eigenbasis_of(HHe):

    operators = []
    rates = []
    for trt in transfer_times:
        ff = trt[0][0]  
        ii = trt[0][1]
        rt = 1.0/trt[1]
        if ii > ff:
            operators.append(qr.qm.ProjectionOperator(ff, ii, dim=e1_dim))
            rates.append(rt)
        else:
            raise Exception("Only downhill rates should"+
                            " be specified explicitely")
        # thermal factor for the rates
        DeltE = HHe.data[ii,ii]-HHe.data[ff,ff] 
        expDEkbT = numpy.exp(-DeltE/(kb_intK*77.0))
        rtU = rt*expDEkbT        
        print("ET transfer time (",ff,"<--",ii,"):", trt[1],"fs")
        print("Uphill time      (",ii,"<--",ff,"):", 1.0/rtU,"fs")
        print(" energy gap : ", qr.convert(DeltE,"int","1/cm"), "1/cm")
        if (1.0/rtU) > ignore_slower_than:
            print(" - Ignoring uphill energy transfer slower than",
                  ignore_slower_than/1000," ps.")
        else:
            operators.append(qr.qm.ProjectionOperator(ii, ff, dim=e1_dim))
            rates.append(rtU)
            print("Setting uphill transfer rate")
        

#
# Electronic only system-bath interaction operator
#
sbi_el = qr.qm.SystemBathInteraction(sys_operators=operators, rates=rates)
sbi_el.set_system(agg_el)

#
# System-bath interaction including vibrational states
#
sbi = qr.qm.SystemBathInteraction(sys_operators=operators,
                                  rates=rates)
sbi.set_system(agg)


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


#
# Laboratory setup
#

from quantarhei import LabSetup
from quantarhei.utils.vectors import X #, Y, Z

lab = LabSetup()
lab.set_polarizations(pulse_polarizations=[X,X,X], detection_polarization=X)

#
# We define time axis for propagation
# 
time_axis = qr.TimeAxis(0.0, propagation_N_steps, propagation_dt)

agg_el.diagonalize()

#
#  Absorption Spectrum calculated by the pathway method
#

# initial conditions before excitation
rho0 = agg_el.get_DensityMatrix(condition_type="thermal", temperature=0.0)
# electronic Hamiltonian
ham = agg_el.get_Hamiltonian()
# first order Liouville pathways
pthways = agg_el.liouville_pathways_1(lab=lab, ham=ham, etol=1.0e-5,
                                       verbose=0) 
# absorption spectrum calculator
mac = qr.MockAbsSpectrumCalculator(time_axis, system=agg_el)
mac.bootstrap(rwa=qr.convert(12200.0,"1/cm","int"), 
              shape="Gaussian")
mac.set_pathways(pthways)

# here we calculate the spectrum
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
    print("\nInitial density matrices: "+
          "\nPure electronic vs. vibrational but traced  over vibrations")
    sig0 = agg.trace_over_vibrations(rho0)
    
    #
    # Impulsive excitation with purely electronic system
    #
    rho0e = agg_el.get_DensityMatrix(condition_type="impulsive_excitation")
    rho0e.normalize2()

    table_data = []
    table_data.append(["Traced","Electronic"])

    with qr.eigenbasis_of(HHe):
        for i in range(1, e1_dim):
            traced = "{0:.8f}".format(numpy.real(sig0.data[i,i]))
            electr = "{0:.8f}".format(numpy.real(rho0e.data[i,i]))
            table_data.append([traced, electr])

        table = AsciiTable(table_data)
        print(table.table)            
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


#
# THIS TAKES FEW SECONDS
#

print("\nPropagating the system (electronic only and vibrational)")

#
# Propagation of the whole system including vibrations 
#
print("Including vibrations ...")
t1 = time.time()
rhot = prop.propagate(rho0)
t2 = time.time()
print("... done. Propagation took:", t2-t1, "sec")

#
# Propagation of the electronic only system
#
print("Electronic only ...")
t1 = time.time()
#with qr.eigenbasis_of(HHe):
rhoet = prope.propagate(rho0e)
t2 = time.time()
print("... done. Propagation took:", t2-t1, "sec")


#
# Time dependent density matrix with vibrations can be traced over the vibrations
#
sigt = agg.trace_over_vibrations(rhot)

# we check the shape of its data
#print(rhot.data.shape)


#
# Here we can compare electronic eigenstate populations at different times
#
print("Checking populations at one point:")
Nt = int(propagation_N_steps/2)
tm = time_axis.data[Nt]
print("Populations at time = ", tm)
table_data = []
table_data.append(["Traced","Electronic"])
with qr.eigenbasis_of(HHe):

    for i in range(1, e1_dim):
        traced = "{0:.8f}".format(numpy.real(sigt.data[Nt,i,i]))
        electr = "{0:.8f}".format(numpy.real(rhoet.data[Nt,i,i]))
        table_data.append([traced, electr])

    table = AsciiTable(table_data)
    print(table.table)            


#
# Plot of the dynamics
#
if show_kinetics:
    with qr.eigenbasis_of(HHe):
        plt.figure(1)
        plt.plot(rhot.TimeAxis.data, 
                 numpy.real(sigt.data[:,comparison_elems[0],
                                      comparison_elems[1]]), "-r")
        plt.plot(rhot.TimeAxis.data,
                 numpy.real(rhoet.data[:,comparison_elems[0],
                                       comparison_elems[1]]), "-b")
        #plt.show()
        plt.savefig(os.path.join(pre_out,"fig01_kinetics_electronic_comp.png"))
        plt.close()
    
    with qr.eigenbasis_of(HH):
        plt.figure(1)
        plt.plot(rhot.TimeAxis.data,
                 numpy.real(rhot.data[:,detail_elems[0],detail_elems[1]]))
        #plt.show()
        plt.savefig(os.path.join(pre_out,"fig02_kinetics_all_detail.png"))
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


trD = agg_el.get_TransitionDipoleMoment()


# ## Liouville pathways

#
# Reloading state for restart
#
if restart:
    restart_now = True
else:
    restart_now = False

if restart_now:
    
    print("\nReloading last state ")
    try:
        state = qr.load_parcel("A_saved_state.qrp")
        N_T2 = state[0]  # here we take the saved value
        agg2 = state[3]
        print("... reload succesful")
    except:
        N_T2 = 0
        #restart = False
        restart_now = False
        print("... reload failed; starting from the beginning")

else:
    N_T2 = 0

#
# THIS TAKES FEW TENS OF SECONDS
#
# In future version of Quantarhei, this call will not be needed 
# (it will be done silently elsewhere, when needed)
print("\nDiagonalization of the aggregate representation:")
t1 = time.time()
agg2.diagonalize()
t2 = time.time()
print("Diagonalized in ", t2-t1, "s")


if pure_deph:
    print("Calculating electronic pure dephasing superoperator")
    t1 = time.time()
    p_deph = qr.qm.ElectronicPureDephasing(agg, dtype="Gaussian")
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
    
    if restart:
        # save the eUt for restart
        qr.save_parcel(eUt, os.path.join(pre_out,"eUt.qrp"))
        
    if stop_after_propagation:
        exit()
     
    # we plot selected evolution superoperato elements
    if evol_super_op_elems2plot is not None:    
        for elem in evol_super_op_elems2plot:
            eUt.plot_element(elem)
   
        plt.savefig(os.path.join(pre_out,"element.png"))
    
    print("Finished in ", t2-t1, "sec")
else:
    print("Dynamics will be calculated on fly")


print("Calculating 2D spectra:")

#
# Laser pulse weighting is on the way
#

#spectrum = qr.DFunction(time, values)
# lab.set_pulse_spectra((0,1,2,3), stype="power", spectrum )

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
    
    #
    #  Select t2 time
    #
    T2 = time_so.data[N_T2]
    print("\n***** Calculating spectrum at t2 = ", T2, " fs *****")
        
          
    if not restart_now:

        #
        # Generation of Liouville pathways
        #
    
        rho0 = agg2.get_DensityMatrix(condition_type="thermal",
                                      temperature=0.0)
    
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
        om_low = qr.convert(frequency_interval[0], "1/cm", "int")
        om_up = qr.convert(frequency_interval[1], "1/cm", "int")
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
    
        print("Saving all pathways at t2 = ", T2, " fs")
        prcl = qr.save_parcel(pw,
                              os.path.join(pre_out,
                                           "pathways"+"_"+str(T2)+".qrp"),
                              comment="pathways at t2 = "+str(T2))
        print("... done")
       
        #
        # Save stuff for a restart
        #
        if restart:
            state = [N_T2, cont, eUt, agg2] #, msc]
            prcl = qr.Parcel()
            prcl.set_content(state)
            prcl.save("A_saved_state.qrp")
            
    else:
        
        cont = state[1]
        eUt = state[2]
        eUt.time = time_so 
        #msc = state[4]
        if N_T2+1 < time_so.length:
            print("Setting next step as ", N_T2+1, " i.e. ",
                  time_so.data[N_T2+1], "fs")
            msc.set_next(N_T2+1) 
        
        else:
            print(" ... already calculated")
        cont.t2axis = time_so
        cont.use_indexing_type(time_so)
        restart_now = False
        

    #
    # propagation of the evolution superoperator
    #
    if (eUt_mode == "jit") and (N_T2+1 < time_so.length) :
        print("Propagating dynamics from ", T2,
              " to ", time_so.data[N_T2+1], " ...")
        
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

    
N_T2_pul = int(N_T2/2)
print("Saving 2D spectrum at",
      time_so.data[N_T2_pul], "fs (as an example)")
twod = cont.get_spectrum(time_so.data[N_T2_pul])

#
# When this cell is left, we either retrieve one spectrum from container,
# or calculate a new one from available pathways
#


#
# Plotting an example spectrum
#
    
# plot_window = [11000,13500,11000,13500]

print("Plotting example 2D spectrum at ", time_so.data[N_T2_pul], " fs")
ex2Dfile = "twod_example_at_T2.png"
with qr.energy_units("1/cm"):
    #plt.figure(1)
    twod.plot(window=plot_window, Npos_contours=10,              
              stype="total", spart="real")
    plt.savefig(os.path.join(pre_out, ex2Dfile))
    #plt.close()
print("... saved into: ", ex2Dfile)

#
# If the plotting window is reasonable, we should cut the unnecessary data
# by trimming the 2D spectrum
#
with qr.energy_units("1/cm"): 
    cont.trimall_to(window=plot_window)


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
print("\nCalculating FFT of the 2D maps")
fcont = cont.fft(window=window, dtype="REPH", offset=t_offset)

#
# Have a look which frequencies we actually have
#
Ndat = len(fcont.axis.data)
print("\nNumber of frequency points:", Ndat)
print("In 1/cm they are:")
with qr.energy_units("1/cm"):
    for k_i in range(Ndat):
        print(k_i, fcont.axis.data[k_i])


#
# Current storage of the specta in the container works on string 
# representation of the frequency. As a consequence the container does not
# recognize units of frequency. We print the frequency 
# in 1/cm, but we have to specify it in internal units
#

with qr.frequency_units("1/cm"):
    sp, show_Npoint = fcont.get_nearest(show_omega)

units = "1/cm"
with qr.energy_units(units):
    print("\nPlotting spectrum at frequency:", fcont.axis.data[show_Npoint], units)
    sp.plot(window=plot_window, Npos_contours=30, 
              stype="total", spart="abs")
    fftfile = "twod_fft_map.png"
    sp.savefig(os.path.join(pre_out, fftfile))
    print("... saved into: ", fftfile)

print("\nSaving spectral containers - both time and fft")
spcontfile_fft = "spectra_container_fft.qrp"
spcontfile_tm = "spectra_container_time.qrp"
qr.save_parcel(fcont,os.path.join(pre_out, spcontfile_fft))
qr.save_parcel(cont,os.path.join(pre_out, spcontfile_tm))
print("... saved into: ", spcontfile_fft)
print("           and: ", spcontfile_tm)

# # Check individual pathways here

skip = True

if not skip:
    with qr.energy_units("1/cm"):
        for pathw in pw:
            print(pathw)


    print("The most intensive pathway (as an example)")
    with qr.energy_units("1/cm"):
        pathw = pw[0]
        print(pathw)


print("")
print("######################################################################")
print("#                                                                    #")
print("#        Simulation finished successfully                            #")
print("#                                                                    #")
print("######################################################################")
#
# EOF
#