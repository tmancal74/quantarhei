# -*- coding: utf-8 -*-
"""

    Test of equilibriation by standard Foerster equations
    in secular approximation. This should lead to cannonical equilibrium
    


"""
from aloe import step
from aloe import world

import numpy

from quantarhei.dev.feature import FeatureFileGenerator
from quantarhei.dev.feature import match_number 

from quantarhei import Molecule
from quantarhei import Aggregate
from quantarhei import CorrelationFunction
from quantarhei import TimeAxis
from quantarhei import ReducedDensityMatrix
from quantarhei import energy_units
from quantarhei import eigenbasis_of

feature = """
#
#  Some types of relaxation theories result in canonical equilibrum
#
#
#

Feature: Some relaxation theories result in cannonical equilibrium 

    Theories such as Redfield relaxation tensor etc. should result in
    cannonical equilibrium

Scenario Outline: Foerster rates for a small chain of molecules

"""

example_1 = """
    Examples:
       | temp | time_step | nsteps | matsu  | atol   |
       | 300  | 1.0       | 10000  | 100    |  0.01  |
       | 200  | 1.0       | 10000  | 100    |  0.02  |
       | 100  | 1.0       | 10000  | 100    |  0.02  |
       | 50   | 1.0       | 10000  | 100    |  0.02  |

"""
