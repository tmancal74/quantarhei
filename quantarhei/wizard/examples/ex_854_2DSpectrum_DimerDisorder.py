# -*- coding: utf-8 -*-

#
#
#  This example assumes 2D spectra of a dimer with a range of energy gaps
#  to be precalculated in a designated directoryß
#
#
#
import os
import numpy
import quantarhei as qr


dname = "sim_up"


#
# Transition energy for which the spectra are calculated
#
Ecalc = 10000.0


#
# load all spectra
#
#for de in des:
#    fname = op.path.join(dname,"twod_E0=1000_dE="+str(de)+".qrp")
    
cont = qr.load_parcel(os.path.join(dname, "cont_p_re.qrp"))
   
spects = []
for ii in range(cont.length()):
    sp = cont.get_spectrum(ii)
    prms = sp.get_log_params()
    dE = prms["dE"]
    
    print("dE =", dE)

de_min = cont.get_spectrum(0).get_log_params()["dE"]
de2 = cont.get_spectrum(1).get_log_params()["dE"]
de_step = de2 - de_min
de_N = cont.length()

print(de_min)
print(de_step)
print(de_N)

#
# These are the values of energy gap that the spectra represent
#
des = qr.ValueAxis(de_min, de_N, de_step)
#print(des.data)


#
# parameters of the disorder
#
width_dis = numpy.zeros(2, dtype=qr.REAL)
width_dis[0] = 100.0 # cm-1
width_dis[1] = 100.0 # cm-1

#
# Average values of the monomer energies
#
E1 = 12500.0
E2 = 12300.0


def weight(E, E0,  width):
    """Weighting function characterizing the disorder
    
    
    """
    return numpy.exp(-numpy.sqrt(numpy.log(2.0))*((E-E0)/width)**2)


#
# Integration parameters
#
de = de_step
N1 = 3
Emin1 = E1 - 490.0
Emin2 = E2 - 490.0

es1 = [Emin1 + de*i for i in range(N1)]
es2 = [Emin2 + de*i for i in range(N1)]

dE_max = es1[N1-1] - es2[0]
dE_min = es2[N1-1] - es1[0]


if (dE_max > des.max) or (dE_min < des.min):
    print("dE max required:", dE_max, " --- dE max available:", des.max)
    print("dE min required:", dE_min, " --- dE min available:",des.min)
    raise Exception("Precalculated spectra are not sufficient to integrated the disorder")

#
# Auxiliary data for storage
#
twod_a = spects[0]
data = numpy.zeros(twod_a.data.shape, dtype=twod_a.data.dtype)

for e1 in es1:
    for e2 in es2:
        
        e_shift = e1 - Ecalc
        DE = e2 - e1
        
        # get the spectrum with DE and shift it by e_shift
        (N_DE, err) = des.locate(DE)
        twod_a = spects[N_DE].deepcopy()
        twod_a.shift_energy(e_shift)
        
        data[:,:] += weight(e1, E1, width_dis[0])* \
                     weight(e2, E2, width_dis[1])*twod_a.data[:,:]
        
        
twod_a.data[:,:] = data[:,:]/(N1**2)

twod_a.plot()




        
    

    

