# -*- coding: utf-8 -*-
import numpy

from ..core.managers import UnitsManaged
from .opensystem import OpenSystem
from ..core.saveable import Saveable
from .. import REAL

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
        
        self.modes = modes
        self.Nmodes = len(self.modes)

        self.sbi = None
        self._has_sbi = False

        self.HH = None    
    
        self._mode_couping_init = False
    

    def set_mode_coupling(self, N, M, val):

        if not self._mode_couping_init:
            self.coupling = numpy.zeros((self.Nmodes, self.Nmodes), dtype=REAL)
            self._mode_couping_init = True

        self.coupling[N,M] = val 
        self.coupling[M,N] = val

        
    
    def build(self):
        from .. import SystemBathInteraction




        ops = []
        rts = []
        mcount = 0
        for md in self.modes:
            #print("Mode:", mcount)
            md.build()


            #
            # Hamiltonian
            #
            Ham = md.get_Hamiltonian()

            # here we have to stretch the Hamiltonian to a new system basis

            #
            # Transition dipole moment
            #
            TrDip = md.get_TransitionDipoleMoment()

            # here we have to stretch the transition dipole to a new system basis


            #
            # Relaxation
            #
            sbi = md.get_SystemBathInteraction()

            nops = sbi.KK.shape[0]
            #print("Number of operators:", nops)
            for ii in range(nops):

                oper = sbi.KK[ii,:,:]

                # here we have to stretch the operator to a new system basis

                ops.append(oper)
                rts.append(sbi.rates[ii])
            
            mcount += 1
    
        sbi = SystemBathInteraction(sys_operators=ops, rates=rts)

        self.sbi = sbi
        self._has_sbi = True
        
