# -*- coding: utf-8 -*-
import numpy

from ..core.managers import UnitsManaged
from .opensystem import OpenSystem
from ..core.saveable import Saveable

class VibrationalSystem(UnitsManaged, Saveable, OpenSystem):
    """ Represents a set of coupled oscillators (possibly unharmonic)


    This class forms the basis of the IR spectroscopy treatment
    in Quantarhei

    Parameters
    ----------

    name : str
        Specifies the name of the system

    modes : list or tuple
        List of modes out of which the systems is built

    """

    def __init__(self, modes=None, name=""):
        
        pass
