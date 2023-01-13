# -*- coding: utf-8 -*-

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

m = qr.Manager()
m.warn_about_basis_change = False 
        
#<<<<<<< HEAD
#sb_reference = qr.BasisReferenceOperator(ham.dim,
#                                      name="site basis reference")


#ham.protect_basis()
with qr.eigenbasis_of(ham):
    
    #RRT = qr.qm.RedfieldRelaxationTensor(ham, sbi, name="Tensor 1")
    #
    #print("\nRelaxation times from the full relaxation tensor")
    #for i in range(1, ham.dim):
    #    for j in range(1, ham.dim):
    #        print(i, "<-", j, ":", 1.0/numpy.real(RRT.data[i,i,j,j]))
        
    print("\nCalculating relaxation rates")
    
    #try: 
    RRM = qr.qm.ModifiedRedfieldRateMatrix(ham, sbi, time)
    #print('The modified redfield rate matrix is:')
    #print(RRM.rates)
    numpy.save('mod_red_rates_pentamer-1_env', RRM.rates)

#except:
    #    pass
    
    #print("\nComparison of the results: ratio of rates")
    #for i in range(1, ham.dim):
    #    for j in range(1, ham.dim):
    #        print(i, "<-", j, ":", RRM.data[i,j]/numpy.real(RRT.data[i,i,j,j]))

    #TDRRM = qr.qm.TDRedfieldRateMatrix(ham, sbi)
    #print("\nRelaxation times from the rate matrix")
    #for i in range(1,ham.dim):
    #    for j in range(1, ham.dim):
    #        print(i, "<-", j, ":", 1.0/TDRRM.data[time.length-1,i,j])

#ham.unprotect_basis()
#=======
#RRM = qr.qm.ModifiedRedfieldRateMatrix(ham, sbi)
print("\nRelaxation times from the rate matrix")

for i in range(0,ham.dim - 1):
    for j in range(0, ham.dim - 1):
         print(i, "<-", j, ":", 1.0/RRM.rates[i,j])

#>>>>>>> upstream/master

if False:
    with qr.eigenbasis_of(ham):
        
        #
        # Evolution of reduced density matrix
        #
    
        prop = qr.ReducedDensityMatrixPropagator(time, ham, RRT)
    
        rho_i = qr.ReducedDensityMatrix(dim=ham.dim, name="Initial DM")
        rho_i.data[3,3] = 1.0
       
        # FIXME: unprotecting does not work correctly
        #RRT.unprotect_basis()
        
        with qr.eigenbasis_of(sb_reference):
            print(" Relaxation time site basis: ", 1.0/RRT.data[1,1,2,2])
    
        RRT.secularize()
        print(" Relaxation time exciton basis: ", 1.0/RRT.data[1,1,2,2])
        rho_t = prop.propagate(rho_i, name="Redfield evolution")
    
        if _show_plots_:
            rho_t.plot(coherences=False)
        
        rho_i1 = qr.ReducedDensityMatrix(dim=ham.dim, name="Initial DM")
        rho_i1.data[3,3] = 1.0   
    
 
    #rho_t.plot(coherences=False)
    
    
    #
    #  Evolution of populations
    #
    
    prop = qr.PopulationPropagator(time, RRM)
    p0 = [i for i in range(ham.dim)]
    p0[3] = 1.0
    pop_t = prop.propagate(p0)
    
    if _show_plots_:
        import matplotlib.pyplot as plt
        plt.plot(time.data, pop_t[:,3],'--r')
        plt.show()
    
#    #print(RRM.data[2,3])
#    #with eigenbasis_of(ham):
#    #    print(RRT.data[2,2,3,3])
    
