# -*- coding: utf-8 -*-

import unittest
import numpy
import tempfile

"""
*******************************************************************************


    Tests of the quantarhei.Molecule class


*******************************************************************************
"""

from quantarhei import Molecule
from quantarhei import Mode
from quantarhei import CorrelationFunction
from quantarhei import TimeAxis
from quantarhei import eigenbasis_of, energy_units


from quantarhei.core.units import kB_intK
from quantarhei import REAL

import quantarhei as qr

class TestOpenSystem_Molecule(unittest.TestCase):

    
    def setUp(self):
        
        import quantarhei as qr
        
        Ng = 1  # number of vibrational states per mode in the electronic ground state
        Ne = 3  # number of vibrational states per mode in the first electronic excited state
        Nf = 3  # number of vibrational states per mode in the second electronic excited state
        
        # use or not the second mode
        second_mode = False
        
        with qr.energy_units("1/cm"):
        
            # three state molecule
            m1 = qr.Molecule([0.0, 12000.0, 14500.0])
        
            # transition dipole moment
            m1.set_dipole((0,1), [0.0, 3.0, 0.0])  # Qy transition
            m1.set_dipole((0,2), [1.5, 0.0, 0.0])  # Qx transition
        
            # first vibrational mode
            mod1 = qr.Mode(1200.0)
            m1.add_Mode(mod1)
            mod1.set_nmax(0,Ng)  # set number of states in the electronic ground state
            mod1.set_nmax(1,Ne)  #     state 1
            mod1.set_nmax(2,Nf)  #     state 2
        
            mod1.set_HR(1,0.1)   # Huang-Rhys factor of the mode in state 1
            mod1.set_HR(2,0.2)  #.    state 2
        
            # second mode is optional
            if second_mode:
        
                mod2 = qr.Mode(1200.0)
                m1.add_Mode(mod2)
                mod2.set_nmax(0,Ng)
                mod2.set_nmax(1,Ne)
                mod2.set_nmax(2,Nf)
                mod2.set_HR(1,0.1)
                mod2.set_HR(2,0.2)

                
            if second_mode:
                #  alpha*Q_1 - beta*Q_2
                alpha = 400.0
                beta = 400.0
                m1.set_diabatic_coupling((1, 2), [alpha, [1,0]])
                m1.set_diabatic_coupling((1, 2), [-beta, [0,1]])
            else:
                # alpha*Q_1
                alpha = 800.0
                m1.set_diabatic_coupling((1,2), [alpha, [1]])

        cfce_params1 = dict(ftype="OverdampedBrownian",
                           reorg=30.0,
                           cortime=50.0,
                           T=300,matsubara=100)
        
        ta = qr.TimeAxis(0.0, 1000, 1.0)
        
        with qr.energy_units("1/cm"):
            cfce = qr.CorrelationFunction(ta, cfce_params1)
        
        m1.set_transition_environment((0,1), cfce)
        m1.set_transition_environment((0,2), cfce)

        
        self.m1 = m1


    def test_initial_condition(self):
        """Initial condition set by excitation"""

        m1 = self.m1
        
        HH = m1.get_Hamiltonian()
        
        with qr.eigenbasis_of(HH):
            #rhoi = m1.get_excited_density_matrix(condition="delta")
            om1 = 10000
            dom = 10.0
            No = 1000
            dat = numpy.zeros(No, dtype=qr.REAL)
            with qr.energy_units("1/cm"):
                ax = qr.FrequencyAxis(om1, 1000, dom)
        
                for io in range(No):
                    om = ax.data[io]
                    dat[io] = numpy.exp(-((om-16000.0)/1000.0)**2)
                    #print(om, dat[io])
        
                # normalize the pulse
                ssum = numpy.sum(dat)*dom
                dat = dat/ssum
                mx = numpy.max(dat)
                dat = dat/mx
        
                spect = qr.DFunction(ax, dat)
        
            #with qr.energy_units("1/cm"):
                print(spect.at(16000.0))
        
        rhoi = m1.get_excited_density_matrix(condition=["pulse_spectrum",spect])





if __name__ == '__main__':
    unittest.main()
