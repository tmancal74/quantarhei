# -*- coding: utf-8 -*-

import os

import matplotlib.pyplot as plt

import quantarhei as qr

dirn = "c01"
pre_in = os.path.join(dirn,"out")
filename = "pathways_20.0.qrp"
Np = 13

ex2Dfile = "test.png"
plot_window = [11000, 13500, 11000, 13500]

t1axis = qr.TimeAxis(0.0, 1000, 1.0)
t2axis = qr.TimeAxis(0.0, 100, 10.0)
t3axis = qr.TimeAxis(0.0, 1000, 1.0)


pws = qr.load_parcel(os.path.join(pre_in,filename))

print(len(pws))

mscal = qr.MockTwoDSpectrumCalculator(t1axis, t2axis, t3axis)
mscal.bootstrap(rwa=qr.convert(12200,"1/cm","int"))

pw = pws[Np]
mscal.set_pathways([pw])

twod = mscal.calculate()

eUt = qr.load_parcel(os.path.join(pre_in,"eUt.qrp"))
oset = qr.load_parcel(os.path.join(dirn,"A_saved_state.qrp"))

with qr.energy_units("1/cm"):
    print(pw)
    twod.plot(window=plot_window, Npos_contours=10,              
              stype="total", spart="real")
    plt.show()
    
ham = oset[4].get_Hamiltonian()
with qr.eigenbasis_of(ham):
    plt.plot(eUt.time.data, eUt.data[:,7,7,6,6])
    plt.plot(eUt.time.data, eUt.data[:,9,9,6,6])
    
    #plt.plot(eUt.time.data, eUt.data[:,7,9,6,6])
    #plt.plot(eUt.time.data, eUt.data[:,9,9,9,9])
    #plt.plot(eUt.time.data, eUt.data[:,6,6,6,6])


plt.savefig(os.path.join(".", ex2Dfile))
