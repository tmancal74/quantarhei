# -*- coding: utf-8 -*-

import quantarhei as qr
import matplotlib.pyplot as plt
import numpy

ta = qr.TimeAxis(0.0, 1000, 1.0)

"""

    Absorption of a monomeric two-level molecule


"""
cfce_params1 = dict(ftype="OverdampedBrownian",
                   reorg=20.0,
                   cortime=100.0,
                   T=300,matsubara=20)

with qr.energy_units("1/cm"):
    m = qr.Molecule("Molecule",[0.0,10000.0])
    cfce1 = qr.CorrelationFunction(ta,cfce_params1)

    m.set_dipole(0,1,[0.0, 1.0, 0.0])
 
    m.set_egcf((0,1),cfce1)   
    
    a1 = qr.AbsSpect(ta,m)
#FIXME: Set this under the units management
    a1.calculate(rwa=10000*qr.core.units.cm2int)

#FIXME: Make AbsSpect plotable and savable (under units management)
    plt.plot(a1.frequency,a1.data,'-r')

print(len(a1.frequency))
print(a1.axis.length)

a1.plot()


save = False
if save:
    ab = numpy.zeros((len(a1.data),2))
    for kk in range(len(a1.data)):
        ab[kk,0] = a1.frequency[kk]
        ab[kk,1] = a1.data[kk]
    
    numpy.savetxt("absfile",ab)
    
    
"""

    Absorption of a simple trimeric aggregate of two-level molecules


"""