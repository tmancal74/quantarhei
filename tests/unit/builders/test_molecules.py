# -*- coding: utf-8 -*-

import unittest
import numpy
#import h5py
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

class TestMolecule(unittest.TestCase):
    """Basic tests for the Molecule class
    
    
    """
    
    def setUp(self):
        self.en = [0.0, 1.0, 2.0]
        self.m = Molecule(name="Molecule",elenergies=self.en)   
        time = TimeAxis(0,1000,1.0)
        params = dict(ftype="OverdampedBrownian", reorg=20, cortime=100,
                      T=300)
        with energy_units("1/cm"):
            fc = CorrelationFunction(time,params)
        self.m.set_transition_environment((0,1),fc)  
        self.fc = fc
        
            
    
    def test_Molecule_instantiation(self):
        """(Molecule) Testing Molecule instantiation
        
        
        """

        
        self.assertEqual(self.m.name,"Molecule")
        for i in range(2):
            self.assertEqual(self.m.elenergies[i],self.en[i])
            
            
    def test_get_Hamiltonian(self):
        """(Molecule) Testing that Molecule returns correct Hamiltonian 
        
        """
        
        H = self.m.get_Hamiltonian()

        h = numpy.zeros((3,3),dtype=REAL)
        h[1,1] = 1.0
        h[2,2] = 2.0
        
        self.assertTrue(numpy.allclose(H.data,h))

    def test_saving_of_molecule(self):
        """(Molecule) Testing the saving capability of the Molecule class
        
        
        """
        use_temporary_file = True
        
        with energy_units("1/cm"):
            mod = Mode(frequency=150)
            mod1 = Mode(frequency=100)
            
        m2 = Molecule(elenergies=[0.0, 2.0])
        m2.add_Mode(mod)
        m2.add_Mode(mod1)

        
        if use_temporary_file: 
            
            #drv = "core"
            #bcs = False
        
            #with h5py.File('tempfile.hdf5', 
            #               driver=drv, 
            #               backing_store=bcs) as f:
            with tempfile.TemporaryFile() as f:
                                             
                #self.m.save_as(f,"Molecule")
                self.m.save(f, test=True)
                
                # reread it
                m = Molecule()
                #m.load_as(f,"Molecule")
                m = m.load(f, test=True)

        else:

            #with h5py.File('tempfile.hdf5') as f: 
            with open('tempfile.qrp', 'wb') as f:                                          
                #self.m.save_as(f,"Molecules")
                self.m.save(f)
            
            #with h5py.File('tempfile.hdf5') as f:
            with open('tempfile.qrp', 'rb') as f:
                m = Molecule()
                #m.load_as(f,"Molecules")
                m = m.load(f)
            

        self.assertEqual(self.m.name, m.name)
        self.assertEqual(self.m.nel, m.nel)
        numpy.testing.assert_array_equal(self.m.elenergies, m.elenergies)
        
        numpy.testing.assert_array_equal(
                self.m.get_transition_environment((0,1)).data, self.fc.data)
        
        #with h5py.File('tempfile.hdf5', 
        #               driver=drv, 
        #               backing_store=bcs) as f:
        with tempfile.TemporaryFile() as f:
                                         
            #self.m.save_as(f,"Molecule")
            m2.save(f, test=True)
            
            # reread it
            m3 = Molecule()
            #m.load_as(f,"Molecule")
            m3 = m3.load(f, test=True)

        self.assertEqual(m2.name, m3.name)
        self.assertEqual(m2.nel, m3.nel)
        numpy.testing.assert_array_equal(m2.elenergies, m3.elenergies)
        
        self.assertEqual(m2.get_Mode(0).get_energy(0), mod.get_energy(0))
        self.assertEqual(m2.get_Mode(1).get_energy(0), mod1.get_energy(0))
        


class TestMoleculeVibrations(unittest.TestCase):
    
    def setUp(self):
        
        en = [0.0, 1.0]
        self.m = Molecule(name="Molecule",elenergies=en)

        mod1 = Mode(frequency=0.1)
        self.m.add_Mode(mod1)
        mod1.set_nmax(0,3)
        mod1.set_nmax(1,3)
    
        en2 = [0.0,0.1,0.1]
        self.m2 = Molecule(name="AdMolecule",elenergies=en2)    
        self.m2.set_adiabatic_coupling(1,2,0.02)
        mod2 = Mode(frequency=0.01)
        self.m2.add_Mode(mod2)
        mod2.set_nmax(0,3)
        mod2.set_nmax(1,3)
        mod2.set_nmax(2,3)
        

    def test_molecule_with_vibrations_1(self):
        """(Molecule) Testing hamiltonian of a two-level molecule with one mode
        
        """
        
        H1 = self.m.get_Hamiltonian()
        H1expected = numpy.diag([0.05, 0.15, 0.25, 1.05, 1.15, 1.25])        
        self.assertTrue(numpy.allclose(H1expected,H1.data))
        
        mod1 = self.m.get_Mode(0)
        mod1.set_nmax(1,2)
        H2 = self.m.get_Hamiltonian()
        H2expected = numpy.diag([0.05, 0.15, 0.25, 1.05, 1.15])        
        self.assertTrue(numpy.allclose(H2expected,H2.data))  
                
        
        
    def test_thermal_density_matrix_0_temp(self):
        """(Molecule) Thermal density matrix: molecule with one mode, zero temperature
        
        """
        
        rho_eq = self.m.get_thermal_ReducedDensityMatrix()

        # Check if the matrix is diagonal        
        dat = rho_eq.data.copy()
        for i in range(dat.shape[0]):
            dat[i,i] = 0.0
        zer = numpy.zeros(dat.shape,dtype=REAL)
        self.assertTrue(numpy.allclose(dat,zer))
        
        # Check if diagonal is a thermal population
        pop = numpy.zeros(dat.shape[0])
        # get temperature from the molecule
        T = self.m.get_temperature()
        self.assertEqual(T,0.0)        
        
        # get density matrix
        rho = self.m.get_thermal_ReducedDensityMatrix()
        rpop = rho.get_populations()
        
        if numpy.abs(T) < 1.0e-10:
            pop[0] = 1.0
            
        else:
            # thermal populations
            pass

        self.assertTrue(numpy.allclose(pop,rpop))
        
        
    def test_thermal_density_matrix_finite_temp(self):
        """(Molecule) Thermal density matrix: molecule with one mode, finite temperature
        
        """
        
        timeAxis = TimeAxis(0.0,1000,1.0)
        params = {"ftype":"OverdampedBrownian",
                  "T":300,
                  "reorg":30,
                  "cortime":100}
        with energy_units("1/cm"):
            cfcion = CorrelationFunction(timeAxis,params)
            self.m.set_transition_environment((0,1),cfcion)
        
        rho_eq = self.m.get_thermal_ReducedDensityMatrix()

        # Check if the matrix is diagonal        
        dat = rho_eq.data.copy()
        for i in range(dat.shape[0]):
            dat[i,i] = 0.0
        zer = numpy.zeros(dat.shape,dtype=REAL)
        self.assertTrue(numpy.allclose(dat,zer))
        
        # Check if diagonal is a thermal population
        pop = numpy.zeros(dat.shape[0])
        # get temperature from the molecule
        T = self.m.get_temperature()
        self.assertTrue(T==300.0)

        
        # get density matrix
        rpop = rho_eq.get_populations()
        
        # get Hamiltonian
        H = self.m.get_Hamiltonian()       

        if numpy.abs(T) < 1.0e-10:
            pop[0] = 1.0
            
        else:
            # thermal populations
            psum = 0.0
            for n in range(pop.shape[0]):
                pop[n] = numpy.exp(-H.data[n,n]/(kB_intK*T))
                psum += pop[n]
            pop *= 1.0/psum
            
        self.assertTrue(numpy.allclose(pop,rpop))
        
        
    def test_thermal_density_matrix_finite_temp_nondiag(self):
        """(Molecule) Thermal density matrix: finite temperature, non-diagonal Hamiltonian
        
        """
        
        timeAxis = TimeAxis(0.0,1000,1.0)
        params = {"ftype":"OverdampedBrownian",
                  "T":300.0,
                  "reorg":30,
                  "cortime":100}
        cfcion = CorrelationFunction(timeAxis,params)
        self.m2.set_egcf((0,1),cfcion)
        self.m2.set_egcf((0,2),cfcion)
        
        rho_eq = self.m2.get_thermal_ReducedDensityMatrix()

        pop = numpy.zeros(rho_eq._data.shape[0])
        # get temperature from the molecule
        T = self.m2.get_temperature()
        self.assertTrue(numpy.abs(T-300.0)<1.0e-10)
 
        # get Hamiltonian
        H = self.m2.get_Hamiltonian() 
        
        with eigenbasis_of(H):
            # get density matrix

            rpop = rho_eq.get_populations()


            if numpy.abs(T) < 1.0e-10:
                pop[0] = 1.0
            
            else:
                # thermal populations
                psum = 0.0
                for n in range(pop.shape[0]):
                    pop[n] = numpy.exp(-H.data[n,n]/(kB_intK*T))
                    psum += pop[n]
                pop *= 1.0/psum

        
        self.assertTrue(numpy.allclose(pop,rpop))
      

class TestMoleculeMultiVibrations(unittest.TestCase):

    
    def setUp(self):
        
        import quantarhei as qr
        
        Ng = 1  # number of vibrational states per mode in the electronic ground state
        Ne = 3  # number of vibrational states per mode in the first electronic excited state
        Nf = 3  # number of vibrational states per mode in the second electronic excited state
        
        # use or not the second mode
        second_mode = True
        
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

        m2 = m1.deepcopy()

        cfce_params1 = dict(ftype="OverdampedBrownian",
                           reorg=30.0,
                           cortime=50.0,
                           T=300,matsubara=100)
        
        ta = qr.TimeAxis(0.0, 1000, 1.0)
        
        with qr.energy_units("1/cm"):
            cfce = CorrelationFunction(ta, cfce_params1)
        
        with qr.energy_units("1/cm"):
            self.menv = CorrelationFunction(ta, cfce_params1)  
            self.menv2 = CorrelationFunction(ta, cfce_params1) 
        
        
        m1.set_transition_environment((0,1), cfce)
        m1.set_transition_environment((0,2), cfce)

        self.m2 = m2
        self.m1 = m1
        
        
        
        # a very simple molecule to test one mode relaxation
        N3_g = 3
        N3_e = 3
        with qr.energy_units("1/cm"):
            m3 = Molecule([0.0, 10000.0])
            mod3 = Mode(300.0)
            m3.add_Mode(mod3)
            mod3.set_nmax(0,N3_g)
            mod3.set_nmax(1,N3_e)
            mod3.set_HR(1,0.1)
        
        self.N3 = N3_g + N3_e
        self.m3 = m3

        # a very simple molecule to test two mode relaxation
        N3_g = 3
        N3_e = 3
        with qr.energy_units("1/cm"):
            m3 = Molecule([0.0, 10000.0])
            mod3 = Mode(300.0)
            m3.add_Mode(mod3)
            mod3.set_nmax(0,N3_g)
            mod3.set_nmax(1,N3_e)
            mod3.set_HR(1,0.1)
            mod4 = Mode(200.0)
            m3.add_Mode(mod4)
            mod4.set_nmax(0,N3_g)
            mod4.set_nmax(1,N3_e)
            mod4.set_HR(1,0.1) 
            
        self.N4 = N3_g**2 + N3_e**2
        self.m4 = m3
        
        # a tree state molecule to test two mode relaxation
        N3_g = 3
        N3_e = 3
        N3_f = 2
        with qr.energy_units("1/cm"):
            m5 = Molecule([0.0, 10000.0, 11000.0])
            mod3 = Mode(300.0)
            m5.add_Mode(mod3)
            mod3.set_nmax(0,N3_g)
            mod3.set_nmax(1,N3_e)
            mod3.set_nmax(2,N3_f)
            mod3.set_HR(1,0.1)
            mod3.set_HR(2,0.1)
            mod4 = Mode(200.0)
            m5.add_Mode(mod4)
            mod4.set_nmax(0,N3_g)
            mod4.set_nmax(1,N3_e)
            mod4.set_nmax(2,N3_f)
            mod4.set_HR(1,0.1) 
            mod4.set_HR(2,0.1)
            
        self.N5 = N3_g**2 + N3_e**2 + N3_f**2
        self.m5 = m5
        
        # a tree state molecule to test two mode relaxation
        N3_g = 1
        N3_e = 3
        N3_f = 2
        with qr.energy_units("1/cm"):
            m6 = Molecule([0.0, 10000.0, 11000.0])
            mod3 = Mode(800.0)
            m6.add_Mode(mod3)
            mod3.set_nmax(0,N3_g)
            mod3.set_nmax(1,N3_e)
            mod3.set_nmax(2,N3_f)
            mod3.set_HR(1,0.1)
            mod3.set_HR(2,0.2)
            mod4 = Mode(1200.0)
            m6.add_Mode(mod4)
            mod4.set_nmax(0,N3_g)
            mod4.set_nmax(1,N3_e)
            mod4.set_nmax(2,N3_f)
            mod4.set_HR(1,0.1) 
            mod4.set_HR(2,0.2)
        
            alpha = 400.0
            beta = 400.0
            m6.set_diabatic_coupling((1, 2), [alpha, [1,0]])
            m6.set_diabatic_coupling((1, 2), [-beta, [0,1]])

        m6.set_transition_environment((0,1), cfce)
        m6.set_transition_environment((0,2), cfce)
        
        self.m6 = m6


    def test_Hamiltonian_etc(self):
        """(Molecule) Testing multimode Molecule 
        
        """

        m1 = self.m1
        

        HH = m1.get_Hamiltonian()
        dip = m1.get_TransitionDipoleMoment()      
        sbi = m1.get_SystemBathInteraction()
        
        self.assertTrue(True)


    def test_potentials(self):
        """(Molecule) Testing multimode Molecule 
        
        """

        exp1 = numpy.array([ 0.45207638,  2.92411307,  3.59543216], 
                           dtype=REAL)
        exp2 = numpy.array([ 0.45207638,  2.9598406,  3.55970463], 
                           dtype=REAL)
        
        m1 = self.m1

        dq = 0.1
        points = [ii*dq - 3.0 for ii in range(60)]
        pots = m1.get_potential_1D(0, points)
        
        pI = numpy.array(pots[0][10], dtype=REAL)
        p0 = numpy.array(pots[1][10], dtype=REAL)

        #print("Coupled   :", pI)
        #print("Un-coupled:", p0)
        
        numpy.testing.assert_allclose(exp1, pI)
        numpy.testing.assert_allclose(exp2, p0)

        pots = m1.get_potential_2D([1,1], [points, points])
        
        # FIXME: Assertions for 2D potential


    def test_mode_environment(self):
        """(Molecule) Testing mode environment 
        
        """
        
        m2 = self.m2
        
        # correlation function must be specified
        with self.assertRaises(Exception) as context:
            m2.set_mode_environment(0, 0)
            
        self.assertTrue("Correlation function not specified." 
                        in str(context.exception))
        
        # this molecule only has modes 0 and 1, i.e. 2 modes
        with self.assertRaises(Exception) as context:
            
            m2.set_mode_environment(2, corfunc=self.menv)

        # mode environments expected in shape (2,3)
        self.assertTrue(m2._has_mode_env.shape == (2,3))

        # this is how the setting should work
        m2.set_mode_environment(0, 1, corfunc=self.menv)
        
        # we can retrieve the correlation function
        menv = m2.get_mode_environment(0, 1)
        self.assertIs(menv, self.menv)
        



    def test_mode_relaxation_one_no_coupling(self):
        """(Molecule) Testing mode relaxation, zero diabatic coupling 
        
        """
        
        expected = {(0,1):213.157173558, 
                    (1,2):106.578586779,
                    (3,4):159.116760287,
                    (4,5):49.3298076331}
        
        m3 = self.m3
        m3.set_mode_environment(0, 1, corfunc=self.menv)
        m3.set_mode_environment(0, 0, corfunc=self.menv)

        HH = m3.get_Hamiltonian()
        
        self.assertEqual(HH.dim, self.N3)
        
        sbi = m3.get_SystemBathInteraction()
        
        time = sbi.TimeAxis
        RT, ham = m3.get_RelaxationTensor(time,relaxation_theory="stR")
        
        #for ii in range(ham.dim):
        #    for jj in range(ii+1,ham.dim):
        #        if numpy.abs(numpy.real(RT.data[ii,ii,jj,jj])) > 1.0e-10:
        #            print(ii,jj, 1.0/numpy.real(RT.data[ii,ii,jj,jj]))
        
        for key in expected.keys():
            ii = key[0]
            jj = key[1]
            self.assertAlmostEqual(expected[key],
                                   1.0/numpy.real(RT.data[ii,ii,jj,jj]))


    def test_mode_relaxation_two_no_coupling(self):
        """(Molecule) Testing mode relaxation, zero diabatic coupling 
        
        """

        expected = {(0, 3): 213.157173558,
                    (1, 4): 213.157173558,
                    (2, 5): 213.157173558,
                    (3, 6): 106.578586779,
                    (4, 7): 106.578586779,
                    (5, 8): 106.578586779,
                    (9, 12): 159.116760287,
                    (10, 13): 159.116760287,
                    (11, 14): 159.116760287,
                    (12, 15): 49.3298076331,
                    (13, 16): 49.3298076331,
                    (14, 17): 49.3298076331}  
        
        m3 = self.m4
        m3.set_mode_environment(0, 1, corfunc=self.menv)
        m3.set_mode_environment(0, 0, corfunc=self.menv)

        HH = m3.get_Hamiltonian()
        
        self.assertEqual(HH.dim, self.N4)
        
        sbi = m3.get_SystemBathInteraction()
        
        time = sbi.TimeAxis
        RT, ham = m3.get_RelaxationTensor(time,relaxation_theory="stR")
        
        #for ii in range(ham.dim):
        #    for jj in range(ii+1,ham.dim):
        #        if numpy.abs(numpy.real(RT.data[ii,ii,jj,jj])) > 1.0e-10:
        #            print(ii,jj, 1.0/numpy.real(RT.data[ii,ii,jj,jj]))
        
        for key in expected.keys():
            ii = key[0]
            jj = key[1]
            self.assertAlmostEqual(expected[key],
                                   1.0/numpy.real(RT.data[ii,ii,jj,jj]))





    def test_mode_relaxation_two_with_coupling(self):
        """(Molecule) Testing two mode relaxation, non-zero diabatic coupling 
        
        """

        expected = {(0, 1): 131.232945781,
                    (0, 3): 213.157173558,
                    (1, 2): 65.6164728903,
                    (1, 4): 213.157173558,
                    (2, 5): 213.157173558,
                    (3, 4): 131.232945781,
                    (3, 6): 106.578586779,
                    (4, 5): 65.6164728903,
                    (4, 7): 106.578586779,
                    (5, 8): 106.578586779,
                    (6, 7): 131.232945781,
                    (7, 8): 65.6164728903,
                    (9, 10): 113.383896941,
                    (9, 12): 159.116760287,
                    (10, 11): 42.4682638296,
                    (10, 13): 159.116760287,
                    (11, 14): 159.116760287,
                    (12, 13): 113.383896941,
                    (12, 15): 49.3298076331,
                    (13, 14): 42.4682638296,
                    (13, 16): 49.3298076331,
                    (14, 17): 49.3298076331,
                    (15, 16): 113.383896941,
                    (16, 17): 42.4682638296,
                    (18, 19): 90.1537254783,
                    (18, 20): 107.689726039,
                    (19, 21): 107.689726039,
                    (20, 21): 90.1537254783}
       
        m3 = self.m5
        
        m3.set_mode_environment(0, 0, corfunc=self.menv)
        m3.set_mode_environment(0, 1, corfunc=self.menv)
        m3.set_mode_environment(0, 2, corfunc=self.menv)

        m3.set_mode_environment(1, 0, corfunc=self.menv2)
        m3.set_mode_environment(1, 1, corfunc=self.menv2)
        m3.set_mode_environment(1, 2, corfunc=self.menv2)

        HH = m3.get_Hamiltonian()
        
        self.assertEqual(HH.dim, self.N5)
        
        sbi = m3.get_SystemBathInteraction()
        
        time = sbi.TimeAxis
        RT, ham = m3.get_RelaxationTensor(time,relaxation_theory="stR")
        
        #for ii in range(ham.dim):
        #    for jj in range(ii+1,ham.dim):
        #        if numpy.abs(numpy.real(RT.data[ii,ii,jj,jj])) > 1.0e-10:
        #            print(ii,jj, 1.0/numpy.real(RT.data[ii,ii,jj,jj]))
        
        for key in expected.keys():
            ii = key[0]
            jj = key[1]
            self.assertAlmostEqual(expected[key],
                                   1.0/numpy.real(RT.data[ii,ii,jj,jj]))



    def test_mode_relaxation_two_with_coupling_and_transition_env(self):
        """(Molecule) Testing combined two mode relaxation, non-zero diabatic coupling 
        
        """
        
        expstring = """
            1 2 727.645664866
            1 3 1300.41161036
            1 4 2474.51923893
            1 5 9161.77892604
            1 6 51848.5444547
            1 7 62770.0671089
            1 8 26666.9361527
            1 9 687779.302323
            1 10 69253.9271796
            1 11 509238.207123
            1 12 119099.317925
            1 13 42496794.0566
            2 3 914.750041535
            2 4 1526.70840275
            2 5 1078.43320489
            2 6 2143.06568258
            2 7 9028.60190128
            2 8 2565.27593516
            2 9 115320.662808
            2 10 9045.38022805
            2 11 139005.643117
            2 12 77328.6259544
            2 13 13228559.0343
            3 4 490.340998427
            3 5 2811.83460894
            3 6 913.812554187
            3 7 1930.01035984
            3 8 6014.56629185
            3 9 76761.2479553
            3 10 5449.82294567
            3 11 359601.262571
            3 12 77857.0113772
            3 13 7281029.05953
            4 5 337.676107267
            4 6 2777.89919921
            4 7 2746.13285431
            4 8 1865.48084417
            4 9 8405.09148975
            4 10 3054.19511958
            4 11 25327.2909135
            4 12 49397.5978727
            4 13 504390.461369
            5 6 1222.2979835
            5 7 2386.5635631
            5 8 743.318986105
            5 9 2204.69485916
            5 10 11519.6914091
            5 11 28250.3318315
            5 12 7596.8066054
            5 13 920883.606458
            6 7 936.617525216
            6 8 3359.84870285
            6 9 1223.21613167
            6 10 2655.49448677
            6 11 1316.928777
            6 12 4470.69910839
            6 13 1699441.53152
            7 8 563.566950268
            7 9 3292.6695198
            7 10 478.739335017
            7 11 1574.10334706
            7 12 9123.06660703
            7 13 131043.878527
            8 9 538.568721881
            8 10 2625.23430172
            8 11 3821.24277384
            8 12 3788.2346337
            8 13 53299.0795067
            9 10 1773.62388926
            9 11 629.599551003
            9 12 761.554383221
            9 13 1360.51316022
            10 11 301.561416614
            10 12 2082.45868179
            10 13 27716.6326621
            11 12 430.144460517
            11 13 916.561339635
            12 13 763.200670786
        
        """
        expected = {}
        for line in expstring.splitlines():
            rel = line.strip()
            if len(rel) > 0:
                nmbrs = rel.split(" ")
                tpl = (int(nmbrs[0]), int(nmbrs[1]))
                val = float(nmbrs[2])
                expected[tpl] = val
            
        m3 = self.m6
        
        m3.set_mode_environment(0, 0, corfunc=self.menv)
        m3.set_mode_environment(0, 1, corfunc=self.menv2)
        m3.set_mode_environment(0, 2, corfunc=self.menv)

        m3.set_mode_environment(1, 0, corfunc=self.menv2)
        m3.set_mode_environment(1, 1, corfunc=self.menv)
        m3.set_mode_environment(1, 2, corfunc=self.menv2)

        HH = m3.get_Hamiltonian()
        
        self.assertEqual(HH.dim, 14)
        
        sbi = m3.get_SystemBathInteraction()

        time = sbi.TimeAxis
        RT, ham = m3.get_RelaxationTensor(time,relaxation_theory="stR")
        
        
        with eigenbasis_of(ham):
        #    rmax = 0.0
        #    for ii in range(ham.dim):
        #        for jj in range(ii+1,ham.dim):
        #            rate = numpy.abs(numpy.real(RT.data[ii,ii,jj,jj]))
        #            if rate > rmax:
        #                rmax = rate
        #            if rate > 1.0e-10:
        #                print(ii,jj, 1.0/numpy.real(RT.data[ii,ii,jj,jj]))
                        
        #print("Max rate:", rmax, 1.0/rmax)
            for key in expected.keys():
                ii = key[0]
                jj = key[1]
                self.assertAlmostEqual(expected[key],
                                       1.0/numpy.real(RT.data[ii,ii,jj,jj]),
                                       places=4)

            

if __name__ == '__main__':
    unittest.main()

       
                
        