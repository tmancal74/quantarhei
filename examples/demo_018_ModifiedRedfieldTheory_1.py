# -*- coding: utf-8 -*-

#
#
#       IN CONSTRUCTION  !!!
#
#


#
# Demo settings
#
_show_plots_ = False


import numpy

import quantarhei as qr
from quantarhei.models.modelgenerator import ModelGenerator

print("""
*******************************************************************************
*                                                                             *
*                    Modified Redfield Theory Demo                            *
*                                                                             *                  
*******************************************************************************
""")

Nt = 1000
dt = 1.0
time = qr.TimeAxis(0.0, Nt, dt)


mg = ModelGenerator()
agg = mg.get_Aggregate_with_environment(name="pentamer-1_env",
                                        timeaxis=time)


agg.build()
agg.diagonalize()

sbi = agg.get_SystemBathInteraction()
ham = agg.get_Hamiltonian()
ham.set_name("Hamiltonian")

print(">>> Hamiltonian ")
with qr.energy_units("1/cm"):
    print(ham)

print("""
*******************************************************************************

                    Calculating relaxation tensor
                  
*******************************************************************************
""")

#m = qr.Manager()
#m.warn_about_basis_change = False 
        

#ham.protect_basis()
#with qr.eigenbasis_of(ham):
if True:
        
    print("\nCalculating relaxation rates")
    
    RRM = qr.qm.ModifiedRedfieldRateMatrix(ham, sbi, time)
    RRT = qr.qm.ModRedfieldRelaxationTensor(ham, sbi)

    #numpy.save('mod_red_rates_pentamer-1_env', RRM.rates)
#ham.unprotect_basis()

print("\nRelaxation times from the rate matrix")

#with qr.eigenbasis_of(ham):
if True:
    for i in range(1,ham.dim-1):
        for j in range(1,ham.dim-1):
            #if numpy.abs(RRM.rates[i,j]) > 1.0e-10:
            print(i, "<-", j, ":", 1.0/numpy.real(RRM.data[i,j])) 
                      #, 1.0/numpy.real(RRT.data[i,i,j,j]), numpy.real(RRT.data[i,j,j,j]))
            
            #else:
            #    print(i, "<-", j, ": inf")





if _show_plots_:
    with qr.eigenbasis_of(ham):
        
        #
        # Evolution of reduced density matrix
        #
    
        prop = qr.ReducedDensityMatrixPropagator(time, ham, RRT)
    
        rho_i = qr.ReducedDensityMatrix(dim=ham.dim, name="Initial DM")
        rho_i.data[3,3] = 1.0
       
        rho_t = prop.propagate(rho_i, name="ModifiedRedfield evolution")
    
        if _show_plots_:
            rho_t.plot(coherences=False)
        
        rho_i1 = qr.ReducedDensityMatrix(dim=ham.dim, name="Initial DM")
        rho_i1.data[3,3] = 1.0   
    
    
    #
    #  Evolution of populations
    #
    """
    prop = qr.PopulationPropagator(time, RRM)
    p0 = [i for i in range(ham.dim)]
    p0[3] = 1.0
    pop_t = prop.propagate(p0)
    
    if _show_plots_:
        import matplotlib.pyplot as plt
        plt.plot(time.data, pop_t[:,3],'--r')
        plt.show()
    """
    
