# -*- coding: utf-8 -*-
import numpy

from quantarhei import TimeAxis
from quantarhei.qm import RateMatrix

tt = TimeAxis(0.0, 1000, 0.1)

N = 3 # number of sites

# Rate matrix
KK = RateMatrix(data=numpy.zeros((N,N)))

print(KK.N)


