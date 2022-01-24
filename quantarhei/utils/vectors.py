# -*- coding: utf-8 -*-

import numpy
from  .types import check_numpy_array 

X = (1.0,0.0,0.0)
Y = (0.0,1.0,0.0)
Z = (0.0,0.0,1.0)
Dxy = (1.0,1.0,0.0)/numpy.sqrt(2)
Axy = (-1.0,1.0,0.0)/numpy.sqrt(2)

def normalize2(vec,norm=1.0):
    """ Normalizes a vector to a specified size """
    vec = check_numpy_array(vec)
    vel = numpy.sqrt(numpy.dot(vec,vec))
    out = (vec/vel)*norm
    return out

def norm(vec):
    """Returns the vector norm (scalar product with itself) """
    vel = numpy.sqrt(numpy.dot(vec,vec))
    return vel


