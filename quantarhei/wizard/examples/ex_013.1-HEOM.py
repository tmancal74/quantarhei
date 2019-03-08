# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot as plt
import quantarhei as qr
print(
"""
*******************************************************************************
    
  Demonstration of electronic dephasing times in a molecular dimer

  Author: Tomas Mancal
          Faculty of Mathematics and Physics
          Charles University
          Ke Karlovu 5
          CZ-121 16 Prague 2
          Czech Republic
          
  Email: tomas.mancal@mff.cuni.cz
  
*******************************************************************************
  
  This script should use Quantarhei package ver. 0.0.45 and older
  
""")
print("  Current version of Quantarhei:", qr.Manager().version)

###############################################################################
#
# Simulation settings
#
###############################################################################
_show_plots_ = True
# make a series of calculations with different resonance energies
_repeate_ = True

# fit selected decays to assign decay time
_fit_ = True

###############################################################################
#
# Simulation parameters
#
###############################################################################

# values of resonance coupling in cm^-1
if _repeate_:
    Js = [5.0, 10.0, 20.0, 30.0, 40.0 , 50.0, 
          75.0, 100.0, 150.0, 200.0, 250.0, 300.0]
    #Js = [50.0, 100.0, 150.]
    #Js_nonzero = Js[1:]
else:
    
    Js = [5]

# energies
E1 = 10000.0         # in cm^-1
E2 = 10000.0         # in cm^-1

# bath parameters
lamb = 50.0          # in cm^-1
tc = 50.0            # in fs
Temperature = 300.0  # in K

# number of time steps on the axis
N = 10000
# integration time step
dt = 0.1
# Hierarchy depth
depth = 4
    
#
###############################################################################
#
#
#   SIMULATION CODE BELOW
#
#
###############################################################################
#    
opt_cohs = []
ext_cohs = []
ext_pops = []
sit_pops = []

for J in Js:
    
    print(
    """
*******************************************************************************
*
*      J =""", J, "1/cm\n*",
"""
*******************************************************************************
""")     
    ###########################################################################
    #
    #   Model system definition
    #
    ###########################################################################
    
    #   Two molecules
    with qr.energy_units("1/cm"):
        m1 = qr.Molecule([0.0, E1])
        m2 = qr.Molecule([0.0, E2])
    
    #   Aggregate is built from the molecules    
    agg = qr.Aggregate([m1, m2])
    
    #   Couplings between them are set
    with qr.energy_units("1/cm"):
        agg.set_resonance_coupling(0,1,J)
    
    #   Interaction with the bath is set through bath correlation functions
    timea = qr.TimeAxis(0.0, N, dt)
    cpar1 = dict(ftype="OverdampedBrownian-HighTemperature", reorg=lamb,
                cortime=tc, T=Temperature)
    #cpar2 = dict(ftype="OverdampedBrownian-HighTemperature", reorg=50,
    #            cortime=50, T=300)
    
    with qr.energy_units("1/cm"):
        cfce1 = qr.CorrelationFunction(timea, cpar1)
    #    cfce2 = qr.CorrelationFunction(timea, cpar1)
        
    m1.set_transition_environment((0, 1), cfce1)
    m2.set_transition_environment((0, 1), cfce1)
        
    #    Aggregate is built
    agg.build()
    
    
    ###########################################################################
    #
    #    Definition of the hierarchy
    #
    ###########################################################################
    
    #    Hamiltonian and the system-bath interaction operator is needed to 
    #    define the Kubo-Tanimura hierarchy
    ham = agg.get_Hamiltonian()
    sbi = agg.get_SystemBathInteraction() 
    
    #    We define the hierarchy and the corresponding propagator
    Hy = qr.KTHierarchy(ham, sbi, depth)
    kprop = qr.KTHierarchyPropagator(timea, Hy)
    
    #   Initial density matrix
    rhoi = qr.ReducedDensityMatrix(dim=ham.dim)
    rhon = qr.ReducedDensityMatrix(dim=ham.dim)

    # FIXME: look at a truly realistic excitation
    #if True:
    with qr.eigenbasis_of(ham):
        rhoi.data[2,2] = 0.5
        rhoi.data[1,1] = 0.5
        rhoi.data[1,2] = 0.5
        rhoi.data[2,1] = 0.5
        rhoi.data[0,1] = 0.5
        rhoi.data[1,0] = 0.5
        
        rhon.data[:,:] = rhoi.data[:,:]
        rhon.data[1,2] = 0.0
        rhon.data[2,1] = 0.0
    # density matrix propagation
    rhot = kprop.propagate(rhoi)
    #rhot_n = kprop.propagate(rhon)
    
    #
    # Plotting density matrix in exciton basis
    #
    pop_e = numpy.zeros(timea.length, dtype=qr.REAL)
    pop_s = numpy.zeros(timea.length, dtype=qr.REAL)
    with qr.eigenbasis_of(ham):
        # population in exciton basis
        pop_e[:] = numpy.real(rhot.data[:,2,2])
        if _show_plots_:
            plt.plot(timea.data, rhot.data[:,0,0])
            plt.plot(timea.data, rhot.data[:,1,1],"-k", label="State 1")
            plt.plot(timea.data, rhot.data[:,2,2],"-r", label="State 2")
            #plt.plot(timea.data, rhot_n.data[:,1,1],"-k",
            #         label="State 1 (incoherent ini)")
            #plt.plot(timea.data, rhot_n.data[:,2,2],"-r",
            #         label="State 2 (incoherent ini)")   
            plt.plot(timea.data, numpy.real(rhot.data[:,1,2]),"--b", 
                     label="Coherence 1-2, real part")    
            plt.plot(timea.data, numpy.imag(rhot.data[:,1,2]),"--g", 
                     label="Coherence 1-2, imag part")
            
            plt.legend(loc="upper right", bbox_to_anchor=(0.95, 0.9))
            plt.title("Density Matrix: Excition basis")
            plt.show()
    # population decay in site basis
    pop_s[:] = numpy.real(rhot.data[:,2,2])    
    
    #
    # Plotting density matrix in site basis
    #
    if _show_plots_:
        plt.plot(timea.data, rhot.data[:,0,0])
        plt.plot(timea.data, rhot.data[:,1,1],"-k", label="Site 1")
        plt.plot(timea.data, rhot.data[:,2,2],"-r", label="Site 2")
        plt.plot(timea.data, numpy.real(rhot.data[:,1,2]),"--b",
                 label="Coherence 1-2, real part")    
        plt.plot(timea.data, numpy.imag(rhot.data[:,1,2]),"--g",
                 label="Coherence 1-2, imag part")
        
        #plt.plot(timea.data, numpy.real(rhot.data[:,0,1]),"-m")
        plt.legend(loc="upper right", bbox_to_anchor=(0.95, 0.97))
        plt.title("Density Matrix: Site basis")
        plt.show()
    
    #
    # Plotting optical coherence
    #
    if _show_plots_:
        with qr.eigenbasis_of(ham):
            
            plt.plot(timea.data, numpy.real(rhot.data[:,0,1]),"-m",
                     label="Real part")
            plt.plot(timea.data, numpy.imag(rhot.data[:,0,1]),"--m",
                     label="Imaginary part")
            plt.plot(timea.data, numpy.imag(rhot.data[:,0,1]),"-k",
                     label="Absolute value")
            
            plt.title("Optical coherence: excition basis (RWA)")
            plt.show()
        
    #
    # Plotting coherences in exciton basis
    #
    ext_coh = numpy.zeros(timea.length, dtype=qr.COMPLEX)
    opt_coh = numpy.zeros(timea.length, dtype=qr.COMPLEX)
    with qr.eigenbasis_of(ham):
        om = ham.data[2,2] - ham.data[1,1]
        ext_coh[:] = rhot.data[:,1,2] #*numpy.exp(-1j*om*timea.data)
        om = ham.data[1,1] - ham.data[0,0]
        opt_coh[:] = rhot.data[:,0,1]*numpy.exp(-1j*om*timea.data)
        
        if _show_plots_:
            plt.plot(timea.data, numpy.abs(ext_coh),"-b",
                     label="Electronic coherence")    
            plt.plot(timea.data, numpy.abs(opt_coh),"-m",
                     label="Optical coherence")
            plt.legend()
            plt.title("Absolute value of coherences: exciton basis")
            plt.show()

    #
    # Plotting coherences in site basis
    #        
    om = ham.data[2,2] - ham.data[1,1]
    ext_cohS = rhot.data[:,1,2] #*numpy.exp(-1j*om*timea.data)
    om = ham.data[1,1] - ham.data[0,0]
    opt_cohS = rhot.data[:,0,1]*numpy.exp(-1j*om*timea.data)
    
    if _show_plots_:
        plt.plot(timea.data, numpy.abs(ext_cohS),"-b",
                 label="Electronic coherence")    
        plt.plot(timea.data, numpy.abs(opt_cohS),"-m",
                 label="Optical coherence")
        plt.legend()
        plt.title("Absolute value of coherences: site basis")
        plt.show()

    #
    # Exponential fit of population decay
    #
    if _fit_:
        
        if J != 0.0:
            pop_e_f = qr.DFunction(timea, pop_e)
            popt_p_e = pop_e_f.fit_exponential(guess=[0.5, 1.0/100, 0.3])
            pop_s_f = qr.DFunction(timea, pop_s)
            popt_p_s = pop_s_f.fit_exponential(guess=[0.5, 1.0/100, 0.3])
            if _show_plots_:
                # plot the fit
                with qr.eigenbasis_of(ham):
                    
                    plt.plot(timea.data, pop_e,"-k",
                             label="Population of exciton 2")
                    plt.plot(timea.data, 
                             popt_p_e[0]*numpy.exp(-popt_p_e[1]*timea.data)
                             +popt_p_e[2],"--g", label="Fit of exciton 2")
                    #plt.show()
            
            print("Fitted population decay time (exciton basis): ",
                  1.0/popt_p_e[1])
            print("Fitted population decay time (site basis): ",
                  1.0/popt_p_s[1])
        else:
            popt_p_e = [0.0, 10000.0]
            popt_p_s = [0.0, 10000.0]
            
        if J == 0.0:
            rat = 0.0
            balance = 0.0
        else:
            rat = popt_p_e[1]
            kBT = qr.core.units.kB_int*Temperature
            print("kBT = ", qr.convert(kBT,"int","1/cm"),"cm^-1")
            with qr.eigenbasis_of(ham):
                balance = rat*numpy.exp(-(ham.data[2,2]-ham.data[1,1])/kBT)
            print("Uphill decay time from canonical detailed balance: ",
                  1.0/balance)
        
            rate_predicted = 0.5*balance + 0.5*rat
            print("Prediction for the decoherence time"+
                  " based on population transfer:", 1.0/rate_predicted)
    

        #
        # Exponential fitting of opt_coh and ext_coh
        #
        
        # Define Quantarhei DFunction which can fit itself
        opt_coh_f = qr.DFunction(timea, numpy.abs(opt_coh))
        popt_o = opt_coh_f.fit_exponential(guess=[0.5, 1.0/200.0, 0.0])
        #print("Optimal parameters for optical coherence: ", popt_o)
        
        ext_coh_f = qr.DFunction(timea, numpy.abs(ext_coh))
        popt_e = ext_coh_f.fit_exponential(guess=[0.5, 1.0/300.0, 0.0])
        #print("Optimal parameters for optical coherence: ", popt_e)
        
        if _show_plots_:
            # plot the fit
            with qr.eigenbasis_of(ham):
                
                plt.plot(timea.data, numpy.abs(opt_coh),"-m",
                         label="optical coherence")
                plt.plot(timea.data, 
                         popt_o[0]*numpy.exp(-popt_o[1]*timea.data)
                         +popt_o[2],"--b", label="optical coherence fit")
                plt.plot(timea.data, numpy.abs(ext_coh),"-r",
                         label="electronic coherence")
                plt.plot(timea.data, 
                         popt_e[0]*numpy.exp(-popt_e[1]*timea.data)
                         +popt_e[2],"--k", label="electronic coherence fit")
                plt.legend()
                plt.title("Mono-exponential fits of"+
                          " coherence and population decay")
                plt.show()
            if J != 0.0:  
                plt.plot(timea.data, pop_s,"-r",
                             label="population")
                plt.plot(timea.data, 
                         popt_p_s[0]*numpy.exp(-popt_p_s[1]*timea.data)
                         +popt_p_s[2],"--k", label="population fit")
                plt.legend()
                plt.title("Mono-exponential fits of"+
                          " site basis population decay")
                plt.show()
            
        print("Fitted optical coherence decay time (exciton basis): ",
              1.0/popt_o[1])
        print("Fitted electronic coherence decay time (exciton basis): ",
              1.0/popt_e[1])
    
    
        opt_cohs.append(1.0/popt_o[1])
        ext_cohs.append(1.0/popt_e[1])
        ext_pops.append(1.0/popt_p_e[1])
        sit_pops.append(1.0/popt_p_s[1])
    
#
# summary
#
if _repeate_ and _fit_:
    print("\nSummary of the results")
    upper = 1000
    opt_cohs = numpy.array(opt_cohs)
    ext_cohs = numpy.array(ext_cohs)
    ext_pops = numpy.array(ext_pops)
    plt.plot(Js, opt_cohs, "-k", label="Optical coherence")
    plt.plot(Js, ext_cohs, "-r", label="Electronic coherence")
    plt.plot(Js, ext_pops, "-b", label="Excitonic population")
    plt.plot(Js, sit_pops, "--g", label="Site population")
    plt.title("Comparison of decay times")
    plt.xlabel(r'J [cm$^{-1}$]')
    plt.ylabel(r'Decay time [fs]')
    plt.legend()
    plt.show()
    plt.plot(Js, opt_cohs, "-k", label="Optical coherence")
    plt.plot(Js, ext_cohs, "-r", label="Electronic coherence")
    plt.plot(Js, ext_pops, "-b", label="Excitonic population")
    plt.plot(Js, sit_pops, "--g", label="Site population")
    plt.plot([lamb, lamb],[0, upper],"--k", linewidth=1.0)
    plt.axis([0,300,0,upper])
    plt.title("Comparison of decay times")
    plt.xlabel(r'J [cm$^{-1}$]')
    plt.ylabel(r'Decay time [fs]')
    plt.legend()
    plt.savefig("summary_lamb="+str(lamb)+".png")
    plt.show()
