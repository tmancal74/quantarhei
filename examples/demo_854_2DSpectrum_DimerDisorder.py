# -*- coding: utf-8 -*-

"""


  This example assumes 2D spectra of a dimer with a range of energy gaps
  to be precalculated in a designated directory ("sim_up" by default)




"""

import os
import numpy
import quantarhei as qr



INP = qr.Input("ex_854_2DSpectrum_DimerDisorder.yaml") 
               #math_allowed_in=["E1", "E2", "width_dis"])

if INP.input_file_from_results:
    ifile = os.path.join(INP.dname, INP.input_file)
else:
    ifile = INP.input_file
    
INP_pre = qr.Input(ifile, math_allowed_in=INP.with_math_allowed_in)

# directory with the precalculated results
dname = INP.dname #"sim_up"

#
# Transition energy for which the spectra are calculated
#
Ecalc = INP_pre.E_P #10000.0
calculated_width = INP_pre.max_available_fwhm


#
# load all spectra
#
#for de in des:
#    fname = op.path.join(dname,"twod_E0=1000_dE="+str(de)+".qrp")
    
cont = qr.load_parcel(os.path.join(dname, INP.container_file))
   
spects = []
for ii in range(cont.length()):
    sp = cont.get_spectrum(ii)
    prms = sp.get_log_params()
    dE = prms["dE"]
    spects.append(sp)
    #print("dE =", dE)

de_min = cont.get_spectrum(0).get_log_params()["dE"]
de2 = cont.get_spectrum(1).get_log_params()["dE"]
de_step = de2 - de_min
de_N = cont.length()

print("de_min:", de_min)
print("de_step:", de_step)
print("de_N", de_N)

#
# These are the values of energy gap that the spectra represent
#
des = qr.ValueAxis(de_min, de_N, de_step)
#print(des.data)


#
# parameters of the disorder
#
width_dis = numpy.zeros(2, dtype=qr.REAL)
width_dis[0] = INP.width_dis[0] #50.0 # cm-1
width_dis[1] = INP.width_dis[1] #50.0 # cm-1
how_many_fwhm = INP.how_many_fwhm

#
# Average values of the monomer energies
#
if INP.use_default_values:
    E1 = INP_pre.E_P #10000.0
    E2 = INP_pre.E_B # INP_pre.E_P + INP_pre.center #9500.0
    if INP_pre.special_pair["useit"]:
        E1 = INP_pre.special_pair["E_Pplus"]

else:
    E1 = INP.E1
    E2 = INP.E2

print("E1:", E1)
print("E2:", E2)

def weight(E, E0,  width):
    """Weighting function characterizing the disorder
    
    
    """
    return numpy.exp(-numpy.sqrt(numpy.log(2.0))*((E-E0)/width)**2)


#
# Integration parameters
#
de = de_step
N1 = de_N
N1_1 = int(width_dis[0]*N1/calculated_width)
N1_2 = int(width_dis[1]*N1/calculated_width)
Emin1 = E1 - de*(N1_1-1)/4.0 #490.0
Emin2 = E2 - de*(N1_2-1)/4.0 # 

es1 = [Emin1 + de*i for i in range(int((N1_1-1)/2)+1)]
es2 = [Emin2 + de*i for i in range(int((N1_2-1)/2)+1)]

dE_max = es2[int((N1_1-1)/2)] - es1[0]  
dE_min = es2[0] - es1[int((N1_2-1)/2)]

if (dE_max > des.max) or (dE_min < des.min):
    print("dE max required:", dE_max, " --- dE max available:", des.max)
    print("dE min required:", dE_min, " --- dE min available:",des.min)
    raise Exception("Precalculated spectra are not "+
                    "sufficient to integrated the disorder")

print("dE max:", dE_max)
print("dE min:", dE_min)

#raise Exception()

#
# Auxiliary data for storage
#
twod_a = spects[0]
data = numpy.zeros(twod_a.data.shape, dtype=twod_a.data.dtype)

#
# Integration over disorder
#
for e1 in es1:
    for e2 in es2:
        
        e_shift = e1 - Ecalc
        DE = e2 - e1
        
        # get the spectrum with DE and shift it by e_shift
        (N_DE, err) = des.locate(DE)
        #print(DE, N_DE, weight(e1, E1, width_dis[0])*
        #                weight(e2, E2, width_dis[1]))
        twod_a = spects[N_DE].deepcopy()
        twod_a.shift_energy(qr.convert(e_shift,"1/cm","int"))
        
        data[:,:] += weight(e1, E1, width_dis[0])* \
                     weight(e2, E2, width_dis[1])*twod_a.data[:,:]
        
        
twod_a.data[:,:] = data[:,:]/(N1_1*N1_2)

with qr.energy_units("1/cm"):
    twod_a.plot(spart=qr.part_ABS)
    if INP.save_fig:
        twod_a.savefig(INP.fig_file)

if INP.save_spectrum:
    twod_a.save(INP.spectrum_file)




        
    

    

