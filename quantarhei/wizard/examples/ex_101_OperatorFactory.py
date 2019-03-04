# -*- coding: utf-8 -*-
#import numpy
import quantarhei as qr
from quantarhei.models import HarmonicOscillator as HO


print("""
   
     Demonstration of the operators from operator_factory object

""")
of = HO.operator_factory(N=3)

sh = of.shift_operator(1j)

print(sh[:3,:3])