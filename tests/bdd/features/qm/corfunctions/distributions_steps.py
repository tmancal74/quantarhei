# -*- coding: utf-8 -*-

from aloe import step
from aloe import world

from ...stepslib import temperature

#from quantarhei import BoseEinsteinDistribution


@step('I create Bose-Einstein distribution')
def bose_einstein_dist(self):
    wa = world.wa 
    T = world.temp
    
    nn = BoseEinsteinDistribution(wa,T)
    
    