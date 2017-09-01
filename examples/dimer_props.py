# -*- coding: utf-8 -*-


import numpy


def theta(J,dE):
    return numpy.arctan(2*J/dE)/2.0

def en_shift(J, E1, E2):
    return numpy.sqrt((E1-E2)**2 + 4*(J**2))/2.0

Js = [100.0, 200.0, 300.0, 400.0, 500.0]

E1 = 16000.0
E2 = 18600.0


for J in Js:

    # theta
    the = theta(J,E2-E1)
    # energy shift
    esh = en_shift(J,E1,E2)

    print("| ", numpy.sin(the), " | ", (E1+E2)/2.0 + esh - E2 )

#    print(" + energy =", (E1+E2)/2.0 + esh)
#    print(" - energy =", (E1+E2)/2.0 - esh)
#    print("delta e = ", ) 
