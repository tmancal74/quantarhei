# -*- coding: utf-8 -*-
"""

    Simulation script 
    
    This script performs a set of simulations by calling the function "run"
    specified below. Scroll down below the definition of the "run" function
    for the simulation parameters.

    The default version of the script runs a series of three simulations
    (parameters Ns_d = 1 and Ns_u=1)
    


"""
import time
import datetime
import os
import shutil
import gc

import numpy

import quantarhei as qr

from quantarhei.utils.vectors import X 
import quantarhei.functions as func
from quantarhei.core.units import kB_int

print("\n***** Calculation of material for disorder integration (dimer version) *****")

input_file = "ex_853_RC.yaml"
#input_file = {'E0': 10000.0, 'resonance_coupling': 100.0, 'no_g_vib': 2, 'no_e_vib': 2, 'params': {'HR': 0.01, 'omega': 500.0, 'use_vib': True}, 'location_of_vibrations': 'up', 'append_to_dirname': '_center=600_FWHM=100', 'dip1': [1.5, 0.0, 0.0], 'dip2': [-1.0, -1.0, 0.0], 'rate': '1.0/500.0', 'temperature': 77.0, 't2_N_steps': 100, 't2_time_step': 10.0, 'fine_splitting': 10, 't1_N_steps': 100, 't1_time_step': 10.0, 't3_N_steps': 100, 't3_time_step': 10.0, 'feature_width': 100.0, 'trim_maps_to': [9900, 11500, 9000, 11500], 'omega_uncertainty': 200.0, 'tukey_window_r': 0.3, 'center': 600.0, 'step': 2.0, 'max_available_fwhm': 100.0, 'how_many_fwhm': 2, 'make_movie': False, 'show_plots': False, 'save_containers': False, 'detailed_balance': True, 't2_save_pathways': [50.0, 100.0, 200.0, 300.0], 'copy_input_file_to_results': True, '_math_allowed_in': ['E0', 'resonance_coupling', 'rate', ['params', ['HR', 'omega', 'rate']], 'center', 'step', 'max_available_fwhm', 'how_many_fwhm', 't2_save_pathways']}
INP = qr.Input(input_file, show_input=True) #, 
               #math_allowed_in =["E0", 
               #                  ["params", ["HR", "omega", "rate"]] ])

make_movie = INP.make_movie
show_plots = INP.show_plots
save_containers = INP.save_containers
detailed_balance = INP.detailed_balance



def run(omega, HR, dE, JJ, rate, E0, vib_loc="up", use_vib=True,
        stype=qr.signal_REPH, make_movie=False, save_eUt=False,
        t2_save_pathways=[], dname=None, trimer=None):
    """Runs a complete set of simulations for a single set of parameters
    
    
    
    """
    if dname is None:
        dname = "sim_"+vib_loc
        
    use_trimer =  trimer["useit"]
        
    #
    #  FIXED PARAMETERS
    #
    if use_trimer:
        dip1 = INP.special_pair["dip1"]
        dip3 = INP.special_pair["dip2"]
    else:
        dip1 = INP.dip_P # [1.5, 0.0, 0.0]
        
    dip2 = INP.dip_B # [-1.0, -1.0, 0.0]
    width = INP.feature_width # 100.0
    #rate = 1.0/50.0
    
    normalize_maps_to_maximu = False
    trim_maps = False
    
    units = "1/cm"
    with qr.energy_units(units):
        
        data_descr = "_dO="+str(dE-omega)+"_HR="+str(HR)+"_J="+str(JJ)
        
        if use_vib:
            sys_char = "_vib"
        else:
            sys_char = "_ele"
        data_ext = sys_char+".png"
        obj_ext = sys_char+".qrp"
        
    #raise Exception()
    
    # parameters of the SP
    if use_trimer:
        E2 = trimer["E_Pminus"]
        epsa = (E0+E2)/2.0
        DE = trimer["DE"]
        J2 = 0.5*numpy.sqrt(((E0-E2)**2)-(DE**2))
        ESP2 = epsa + DE/2.0
        ESP1 = epsa - DE/2.0
        rate_3 = trimer["rate"]
        
    use_rate_3 = True
        
    #
    #   Model system is a dimer of molecules
    #
    with qr.energy_units("1/cm"):
        if not use_trimer:
            
            mol1 = qr.Molecule([0.0, E0])
            mol2 = qr.Molecule([0.0, E0+dE])
            
            print("Monomer 1 energy:", E0)
            print("Monomer 2 energy:", E0+dE)
        else:
            mol1 = qr.Molecule([0.0, ESP2])
            mol2 = qr.Molecule([0.0, E0+dE])
            print("Monomer 1 energy:", ESP2)
            print("Monomer 2 energy:", E0+dE)            
            mol3 = qr.Molecule([0.0, ESP1])
            print("Monomer 3 energy:", ESP1)
            
            mol3.set_transition_width((0,1), qr.convert(width, "1/cm", "int"))
            mol3.set_dipole(0,1, dip3)

        
        mod1 = qr.Mode(omega)
        mod2 = qr.Mode(omega)
    
    mol1.set_transition_width((0,1), qr.convert(width, "1/cm", "int"))
    mol1.set_dipole(0,1, dip1)
    
    mol2.set_transition_width((0,1), qr.convert(width, "1/cm", "int"))
    mol2.set_dipole(0,1, dip2)
    
    if use_trimer:
        agg = qr.Aggregate([mol1, mol2, mol3])
    else:
        agg = qr.Aggregate([mol1, mol2])
    
    with qr.energy_units("1/cm"):
        agg.set_resonance_coupling(0,1,JJ)
        if use_trimer:
            agg.set_resonance_coupling(0,2,J2)
            print("SP coupling:", J2)
    
    #
    # Electronic only aggregate
    #
    agg_el = agg.deepcopy()
        
    #
    # if nuclear vibrations are to be added, do it here
    #
    if use_vib:
    
        if vib_loc == "down":
            set_vib = [True, False]
        elif vib_loc == "up":
            set_vib = [False, True]
        elif vib_loc == "both":
            set_vib = [True, True]
        else:
            raise Exception("Unknown location of the vibrations")
            
        if set_vib[0]:
            mol1.add_Mode(mod1)
            mod1.set_nmax(0, INP.no_g_vib)
            mod1.set_nmax(1, INP.no_e_vib)
            mod1.set_HR(1, HR)
        
        if set_vib[1]:
            mol2.add_Mode(mod2)
            mod2.set_nmax(0, INP.no_g_vib)
            mod2.set_nmax(1, INP.no_e_vib)
            mod2.set_HR(1, HR)
    
    agg3 = agg.deepcopy()
    
    agg.build(mult=1)
    agg_el.build(mult=1)
    
    HH = agg.get_Hamiltonian()
    He = agg_el.get_Hamiltonian()
    
    #with qr.energy_units("1/cm"):
    #    print(He)
    
    with qr.energy_units("1/cm"):
        with qr.eigenbasis_of(He):
            Ep_l = He.data[1,1]
            Ep_u = He.data[2,2]

    Ep = numpy.zeros((4,2))
    Ep[0,0] = Ep_l
    Ep[0,1] = Ep_l
    Ep[1,0] = Ep_l
    Ep[1,1] = Ep_u
    Ep[2,0] = Ep_u
    Ep[2,1] = Ep_l
    Ep[3,0] = Ep_u
    Ep[3,1] = Ep_u

        
    #
    # Laboratory setup
    #
    lab = qr.LabSetup()
    lab.set_polarizations(pulse_polarizations=[X,X,X], 
                          detection_polarization=X)
    
    t2_N_steps = INP.t2_N_steps
    t2_time_step = INP.t2_time_step
    time2 = qr.TimeAxis(0.0, t2_N_steps, t2_time_step)
    
    cont_p = qr.TwoDResponseContainer(t2axis=time2)
    cont_m = qr.TwoDResponseContainer(t2axis=time2) 
    #
    # spectra will be indexed by the times in the time axis `time2`
    #
    cont_p.use_indexing_type(time2)
    
    #
    # We define two-time axes, which will be FFTed and will define 
    # the omega_1 and omega_3 axes of the 2D spectrum
    #
    t1_N_steps = INP.t1_N_steps
    t1_time_step = INP.t1_time_step
    t3_N_steps = INP.t3_N_steps
    t3_time_step = INP.t3_time_step
    t1axis = qr.TimeAxis(0.0, t1_N_steps, t1_time_step)
    t3axis = qr.TimeAxis(0.0, t3_N_steps, t3_time_step)
    
    #
    # This calculator calculated 2D spectra from the effective width 
    #
    msc = qr.MockTwoDResponseCalculator(t1axis, time2, t3axis)
    with qr.energy_units("1/cm"):
        msc.bootstrap(rwa=E0, shape="Gaussian")

    #
    # System-bath interaction including vibrational states
    #
    operators = []
    rates = []

    with qr.eigenbasis_of(He):
        if (He.data[2,2] < He.data[1,1]): # or (He.data[3,3]>He.data[2,2]):
            Exception("Electronic states not orderred!")
            
        
        if use_trimer:
            # B -> P
            operators.append(qr.qm.ProjectionOperator(2, 3, dim=He.dim))
            rates.append(rate)
            if use_rate_3:
                # P+ -> P-
                operators.append(qr.qm.ProjectionOperator(1, 2, dim=He.dim))
                rates.append(rate_3)
        else:
            # B -> P
            operators.append(qr.qm.ProjectionOperator(1, 2, dim=He.dim))
            rates.append(rate)        
    
    # include detailed balace
    if detailed_balance:
        with qr.eigenbasis_of(He):
            T = INP.temperature #77.0
            if use_trimer:
                Den = (He.data[3,3] - He.data[2,2])/(kB_int*T)
                operators.append(qr.qm.ProjectionOperator(3, 2, dim=He.dim))
                thermal_fac = numpy.exp(-Den)
                rates.append(rate*thermal_fac)   
                if use_rate_3:
                    Den = (He.data[2,2] - He.data[1,1])/(kB_int*T)
                    operators.append(qr.qm.ProjectionOperator(2, 1, dim=He.dim))
                    thermal_fac_3 = numpy.exp(-Den)
                    rates.append(rate_3*thermal_fac_3)
            else:
                Den = (He.data[2,2] - He.data[1,1])/(kB_int*T)
                operators.append(qr.qm.ProjectionOperator(2, 1, dim=He.dim))
                thermal_fac = numpy.exp(-Den)
                rates.append(rate*thermal_fac)
    
    sbi = qr.qm.SystemBathInteraction(sys_operators=operators, rates=rates)
    sbi.set_system(agg)
    
    #
    # Liouville form for relaxation
    #
    LF = qr.qm.ElectronicLindbladForm(HH, sbi, as_operators=True)
    
    #
    # Pure dephasing
    #
    p_deph = qr.qm.ElectronicPureDephasing(agg, dtype="Gaussian")
    
    
    eUt = qr.qm.EvolutionSuperOperator(time2, HH, relt=LF, pdeph=p_deph,
                                       mode="all")
    eUt.set_dense_dt(INP.fine_splitting)
    
    #
    # We calculate evolution superoperator
    #
    eUt.calculate(show_progress=False)
    
    if save_eUt:
        eut_name = os.path.join(dname, "eUt"+
                                    "_omega2="+str(omega)+data_descr+obj_ext)
        eUt.save(eut_name)
    
    #
    # Prepare aggregate with all states (including 2-EX band)
    #
    agg3.build(mult=2)
    agg3.diagonalize()
    
    pways = dict()
    
    olow_cm = omega-INP.omega_uncertainty/2.0
    ohigh_cm = omega+INP.omega_uncertainty/2.0
    olow = qr.convert(olow_cm, "1/cm", "int")
    ohigh = qr.convert(ohigh_cm, "1/cm", "int")
    
    for t2 in time2.data:
        
        # this could save some memory of pathways become too big
        pways = dict()
        

        twod = msc.calculate_one_system(t2, agg3, eUt, lab, pways=pways,
                                        selection=[["omega2",[olow, ohigh]]])

        #print("t2 =", t2)    
        #print("Number of pathways used for omega2 =",omega,":",
        #      len(pways[str(t2)]))

        if t2 in t2_save_pathways:
            pws_name = os.path.join(dname, "pws_t2="+str(t2)+
                                    "_omega2="+str(omega)+data_descr+obj_ext)
            qr.save_parcel(pways[str(t2)], pws_name) 
       
        cont_p.set_spectrum(twod)

        twod = msc.calculate_one_system(t2, agg3, eUt, lab, pways=pways,
                                        selection=[["omega2",[-ohigh, -olow]]])
    
        #print("Number of pathways used for omega2 =",-omega,":",
        #      len(pways[str(t2)]))
        
        if t2 in t2_save_pathways:
            pws_name = os.path.join(dname, "pws_t2="+str(t2)+
                                    "_omega2="+str(-omega)+data_descr+obj_ext)
            qr.save_parcel(pways[str(t2)], pws_name)
        
        cont_m.set_spectrum(twod)
    
    if make_movie:
        with qr.energy_units("1/cm"):
            cont_p.make_movie("mov.mp4")
     
    #
    # let's not save all the pathways
    #
    #fname = os.path.join("sim_"+vib_loc, "pways.qrp")
    #qr.save_parcel(pways, fname)
    
    fname = os.path.join(dname, "aggregate.qrp")
    agg3.save(fname)
        
    #
    # Window function for subsequenty FFT
    #
    window = func.Tukey(time2, r=INP.tukey_window_r, sym=False)
    
    #
    # FFT with the window function
    #
    # Specify REPH, NONR or `total` to get different types of spectra
    #
    print("Calculating FFT of the 2D maps")
    #fcont = cont.fft(window=window, dtype=stype) #, dpart="real", offset=0.0)
    
    fcont_p_re = cont_p.fft(window=window, dtype=qr.signal_REPH)
    fcont_p_nr = cont_p.fft(window=window, dtype=qr.signal_NONR)
    fcont_p_to = cont_p.fft(window=window, dtype=qr.signal_TOTL)
    
    if normalize_maps_to_maximu:
        fcont_p_re.normalize2(dpart=qr.part_ABS)
        fcont_p_nr.normalize2(dpart=qr.part_ABS)
        fcont_p_to.normalize2(dpart=qr.part_ABS)
    
    fcont_m_re = cont_m.fft(window=window, dtype=qr.signal_REPH)
    fcont_m_nr = cont_m.fft(window=window, dtype=qr.signal_NONR)
    fcont_m_to = cont_m.fft(window=window, dtype=qr.signal_TOTL)

    if normalize_maps_to_maximu:   
        fcont_m_re.normalize2(dpart=qr.part_ABS)
        fcont_m_nr.normalize2(dpart=qr.part_ABS)
        fcont_m_to.normalize2(dpart=qr.part_ABS)
   
    if trim_maps:
        twin = INP.trim_maps_to
        with qr.energy_units("1/cm"):
            fcont_p_re.trimall_to(window=twin)
            fcont_p_nr.trimall_to(window=twin)
            fcont_p_to.trimall_to(window=twin)
        
    show_omega = omega
    
    #
    # Have a look which frequencies we actually have
    #
#    Ndat = len(fcont_re.axis.data)
#    print("\nNumber of frequency points:", Ndat)
#    print("In 1/cm they are:")
#    with qr.energy_units("1/cm"):
#        for k_i in range(Ndat):
#            print(k_i, fcont_re.axis.data[k_i])
    
    with qr.frequency_units("1/cm"):
        sp1_p_re, show_Npoint1 = fcont_p_re.get_nearest(show_omega)
        sp2_p_re, show_Npoint2 = fcont_p_re.get_nearest(-show_omega)
        sp1_p_nr, show_Npoint1 = fcont_p_nr.get_nearest(show_omega)
        sp2_p_nr, show_Npoint2 = fcont_p_nr.get_nearest(-show_omega)
        sp1_p_to, show_Npoint1 = fcont_p_to.get_nearest(show_omega)
        sp2_p_to, show_Npoint2 = fcont_p_to.get_nearest(-show_omega)
        sp1_m_re, show_Npoint1 = fcont_m_re.get_nearest(show_omega)
        sp2_m_re, show_Npoint2 = fcont_m_re.get_nearest(-show_omega)
        sp1_m_nr, show_Npoint1 = fcont_m_nr.get_nearest(show_omega)
        sp2_m_nr, show_Npoint2 = fcont_m_nr.get_nearest(-show_omega)
        sp1_m_to, show_Npoint1 = fcont_m_to.get_nearest(show_omega)
        sp2_m_to, show_Npoint2 = fcont_m_to.get_nearest(-show_omega)    
        

    with qr.energy_units(units):


        if show_plots:
        
            print("\nPlotting and saving spectrum at frequency:", 
                  fcont_p_re.axis.data[show_Npoint1], units)
            
            fftf_1 = os.path.join(dname, "twod_fft"+data_descr+
                                   "_stype=REPH"+"_omega="+str(omega)+data_ext)
            sp1_p_re.plot(Npos_contours=10, spart=qr.part_ABS, 
                          label="Rephasing\n $\omega="+str(omega)+
                          "$ cm$^{-1}$", text_loc=[0.05,0.1], 
                          show_states=[Ep_l, Ep_u, Ep_u+numpy.abs(omega)], 
                          show_diagonal="-k")   
            sp1_p_re.savefig(fftf_1)
            print("... saved into: ", fftf_1)
            fftf_2 = os.path.join(dname, "twod_fft"+data_descr+
                                   "_stype=NONR"+"_omega="+str(omega)+data_ext)
            sp1_p_nr.plot(Npos_contours=10, spart=qr.part_ABS, 
                          label="Non-rephasing\n $\omega="+str(omega)+
                          "$ cm$^{-1}$", text_loc=[0.05,0.1],
                          show_states=[Ep_l, Ep_u, Ep_u+numpy.abs(omega)],
                          show_diagonal="-k")   
            sp1_p_nr.savefig(fftf_2)
            print("... saved into: ", fftf_2)
            fftf_3 = os.path.join(dname, "twod_fft"+data_descr+
                                   "_stype=tot"+"_omega="+str(omega)+data_ext)
            sp1_p_to.plot(Npos_contours=10, spart=qr.part_ABS, 
                          label="Total\n $\omega="+str(omega)+
                          "$ cm$^{-1}$", text_loc=[0.05,0.1],
                          show_states=[Ep_l, Ep_u, Ep_u+numpy.abs(omega)],
                          show_diagonal="-k")   
            sp1_p_to.savefig(fftf_3)
            print("... saved into: ", fftf_3)        
        
        #
        # Point evolutions at the expected peak positions
        #

        if show_plots:

            for ii in range(4):        
                points = fcont_p_re.get_point_evolution(Ep[ii,0], Ep[ii,1],
                                                        fcont_p_re.axis)
                points.apply_to_data(numpy.abs)
                if ii >= 3:
                    points.plot(show=True)
                else:
                    points.plot(show=False)
            
        
            print("\nPlotting and saving spectrum at frequency:", 
                  fcont_m_re.axis.data[show_Npoint2], units)
            fftf_4 = os.path.join(dname, "twod_fft"+data_descr+
                                   "_stype=REPH"+"_omega="+str(-omega)+data_ext)
            sp2_m_re.plot(Npos_contours=10, spart=qr.part_ABS,
                          label="Rephasing\n $\omega="+str(-omega)+
                          "$ cm$^{-1}$", text_loc=[0.05,0.1],
                          show_states=[Ep_l, Ep_u, Ep_u+numpy.abs(omega)],
                          show_diagonal="-k")   
            sp2_m_re.savefig(fftf_4)
            print("... saved into: ", fftf_4)
            fftf_5 = os.path.join(dname, "twod_fft"+data_descr+
                                   "_stype=NONR"+"_omega="+str(-omega)+data_ext)
            sp2_m_nr.plot(Npos_contours=10, spart=qr.part_ABS,
                          label="Non-rephasing\n $\omega="+str(-omega)+
                          "$ cm$^{-1}$", text_loc=[0.05,0.1],
                          show_states=[Ep_l, Ep_u, Ep_u+numpy.abs(omega)],
                          show_diagonal="-k")      
            sp2_m_nr.savefig(fftf_5)
            print("... saved into: ", fftf_5)
            fftf_6 = os.path.join(dname, "twod_fft"+data_descr+
                                   "_stype=tot"+"_omega="+str(-omega)+data_ext)
            sp2_m_to.plot(Npos_contours=10, spart=qr.part_ABS,
                          label="Total\n $\omega="+str(-omega)+
                          "$ cm$^{-1}$", text_loc=[0.05,0.1],
                          show_states=[Ep_l, Ep_u, Ep_u+numpy.abs(omega)],
                          show_diagonal="-k")      
            sp2_m_to.savefig(fftf_6)
            print("... saved into: ", fftf_6)

        if show_plots:
            #
            # Point evolutions at the expected peak positions
            #
            for ii in range(4):        
                points = fcont_p_re.get_point_evolution(Ep[ii,0], Ep[ii,1],
                                                        fcont_m_re.axis)
                points.apply_to_data(numpy.abs)
                if ii >= 3:
                    points.plot(show=True)
                else:
                    points.plot(show=False)
                    
            #points.apply_to_data(numpy.abs)
            #points.plot()

    
    # saving containers
#    fname = os.path.join("sim_"+vib_loc,"fcont_re"+data_descr+obj_ext)
#    print("Saving container into: "+fname)
#    fcont_p_re.save(fname)
#    fname = os.path.join("sim_"+vib_loc,"fcont_nr"+data_descr+obj_ext)
#    print("Saving container into: "+fname)
#    fcont_p_nr.save(fname)
#    fname = os.path.join("sim_"+vib_loc,"fcont_to"+data_descr+obj_ext)
#    print("Saving container into: "+fname)
#    fcont_p_to.save(fname)
    
    if save_containers:
        fname = os.path.join(dname,"cont_p"+data_descr+obj_ext)
        print("Saving container into: "+fname)
        cont_p.save(fname)
        fname = os.path.join(dname,"cont_m"+data_descr+obj_ext)
        print("Saving container into: "+fname)
        cont_m.save(fname)
        
    return (sp1_p_re, sp1_p_nr, sp2_m_re, sp2_m_nr)


###############################################################################
###############################################################################
###############################################################################
#
#   PARAMETERS OF THE SIMULATION
#
###############################################################################  


parms1 = [INP.params] #[strings_2_floats(INP.params, keys=["HR", "omega", "rate"])] 
         #[dict(HR=0.01, omega=500.0, rate=1.0/500.0,
         #           use_vib=True)]
# this is a fix to have rate defined separately from other parameters
parms1[0]["rate"] = INP.rate

#
#
#   MODELS (vibrations on different molecules)
#
#
vib_loc = INP.location_of_vibrations
models = [dict(vib_loc=vib_loc)] #, 
#          dict(vib_loc="down"),
#          dict(vib_loc="both")]

#
# t2s at which pathways will be saved
#
t2_save_pathways = INP.t2_save_pathways #[50.0, 100.0, 200.0, 300.0]

trimer = INP.special_pair #INP.trimer

use_trimer =  trimer["useit"]

#
# Here we construct a path through parameters space
#
if use_trimer:
    center = INP.E_B - INP.special_pair["E_Pplus"]
else:
    center = INP.E_B - INP.E_P #600.0

step = INP.step #2.0

max_available_fwhm = INP.max_available_fwhm #100.0
how_many_fwhm = INP.how_many_fwhm #2

#step = 2
Ns_d = int(2.0*how_many_fwhm*max_available_fwhm/step) # 50
Ns_u = int(2.0*how_many_fwhm*max_available_fwhm/step) # 50

vax = qr.ValueAxis(center-Ns_d*step, Ns_d+Ns_u+1, step)


trimer_disorder = False # trimer["disorder"]
#if use_trimer:
#    vax2 = qr.ValueAxis(trimer["center2"]-Ns_d*step, Ns_d+Ns_u+1, step)
    
print("\nSummary of simulation parameters\n")
print("Energy gap values:")
print("Minimal gap =", vax.min)
print("Maximum gap =", vax.max)
#if use_trimer and trimer_disorder:
#    print("Minimal gap (2nd dim) =", vax2.min)
#    print("Maximum gap (2nd dim) =", vax2.max)
print("Number of steps =", vax.length)

#
# Here we specify pairs of parameters (resonance coupling J and energy
# gap \Delta E between the monomers). One could specify an arbitrary 
# "pathway" in the parameters space. Below we specify a line of 
# increasing \Delta E with constant J.
#
ptns = []

single_run = INP.single_realization

if single_run:
    
    ptns.append((INP.resonance_coupling, center, INP.trimer))
    
else:

    if use_trimer:
#        if trimer_disorder:
#            for val in vax.data:
#                for val2 in vax2.data:
#                    ptns.append((INP.resonance_coupling, val, 
#                                 INP.trimer, val2))
#        else:
        for val in vax.data:
            ptns.append((INP.resonance_coupling, val, 
                         INP.special_pair))

    else:
        for val in vax.data:
            ptns.append((INP.resonance_coupling, val, 
                         INP.special_pair))

if use_trimer:
    E0 = INP.special_pair["E_Pplus"]
else:
    E0 = INP.E_P # transition energy (in 1/cm) of the reference monomer


###############################################################################
###############################################################################
###############################################################################

#
# Container for resulting 2D maps
#
cont_p_re = qr.TwoDSpectrumContainer()
cont_p_re.use_indexing_type("integer")
cont_p_nr = qr.TwoDSpectrumContainer()
cont_p_nr.use_indexing_type("integer")
cont_m_re = qr.TwoDSpectrumContainer()
cont_m_re.use_indexing_type("integer")
cont_m_nr = qr.TwoDSpectrumContainer()
cont_m_nr.use_indexing_type("integer")


#
#
#   LOOP OVER SIMULATIONS
#
#

parms = parms1  #parms1+parms2+parms3+parms4+parms5

i_p_re = 0
tags = []

tA = time.time()
at = '{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
print("\nStarting simulation set at", at)
ll = 1
for model in models:
    print("Model no.", ll, "of", len(models))
    vib_loc = model["vib_loc"]
    #stype = model["stype"]
    dname = "sim_"+vib_loc+INP.append_to_dirname
    try:
        os.makedirs(dname)
    except FileExistsError:
        # directory already exists
        pass

    if INP.copy_input_file_to_results:
        if INP._from_file:
            shutil.copy2(input_file, dname)
        else:
            INP.dump_yaml(os.path.join(dname, "input_file.yml"))


    kk = 1
    at = '{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
    print("\nStarting the simulation at:", at)
    Np = len(parms)
    for par in parms:

        print("Run no.", kk, "of", Np)
        omega = par["omega"]
        HR = par["HR"]
        #dE = par["dE"]
        #JJ = par["JJ"]
        use_vib = par["use_vib"]
        rate = par["rate"]

        kp = 1
        Nje = len(ptns)
        for (JJ, dE, trimer) in ptns:
            print("\nCalculating spectra ... (",kp,"of",Nje,") [run ",kk,"of",
                  Np,"]")
            print("JJ =", JJ)
            print("dE =", dE)
            t1 = time.time()
            if Nje == 1:
                save_eUt = True
            else:
                save_eUt = False
            (sp1_p_re, sp1_p_nr, sp2_m_re, sp2_m_nr) = \
            run(omega, HR, dE, JJ, rate, E0, vib_loc, use_vib,
                make_movie=make_movie, save_eUt=save_eUt, 
                t2_save_pathways=t2_save_pathways, dname=dname, trimer=trimer)
            t2 = time.time()
            gc.collect()
            print("... done in",t2-t1,"sec")
        
            params = dict(J=JJ, dE=dE, E0=E0, omega=omega)
            sp1_p_re.log_params(params)
            sp1_p_nr.log_params(params)
            sp2_m_re.log_params(params)
            sp2_m_nr.log_params(params)
            
            cont_p_re.set_spectrum(sp1_p_re, tag=i_p_re)
            cont_p_nr.set_spectrum(sp1_p_nr, tag=i_p_re)
            cont_m_re.set_spectrum(sp2_m_re, tag=i_p_re)
            cont_m_nr.set_spectrum(sp2_m_nr, tag=i_p_re)
            tags.append(i_p_re)
            i_p_re +=1
            kp += 1
            
        kk += 1
    ll += 1
    
tB = time.time()
at = '{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
print("\n... finished simulation set at", at, "in", tB-tA,"sec")

fname = os.path.join(dname, "cont_p_re.qrp")
cont_p_re.save(fname)
fname = os.path.join(dname, "cont_p_nr.qrp")
cont_p_nr.save(fname)
fname = os.path.join(dname, "cont_m_re.qrp")
cont_m_re.save(fname)
fname = os.path.join(dname, "cont_m_nr.qrp")
cont_m_nr.save(fname)


