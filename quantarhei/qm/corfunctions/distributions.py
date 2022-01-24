# -*- coding: utf-8 -*-

import numpy


class BoseEinsteinDistribution:
    
    def __init__(self,freq_axis, temperature):
        kBT = temperature
        self.data = 1.0/(numpy.exp(freq_axis.data/kBT) - 1.0)
        
        