# -*- coding: utf-8 -*-

import quantarhei as qr

from quantarhei.qm.liouvillespace.integrodiff.integrodiff import IntegrodiffPropagator

timea = qr.TimeAxis(0.0, 100, 0.1)
ham = qr.Hamiltonian(data=[[0.0, 1.0], [1.0, 10.0]])

ip = IntegrodiffPropagator(timea, ham)

rhoi= qr.ReducedDensityMatrix(data=[[0.0, 0.0],[0.0, 1.0]])

ip.propagate(rhoi)