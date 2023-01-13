# -*- coding: utf-8 -*-

import unittest
#import h5py
import tempfile
import numpy

"""
*******************************************************************************


    Tests of the quantarhei.Aggregate class


*******************************************************************************
"""

from quantarhei import Aggregate
from quantarhei.builders.aggregate_test import TestAggregate
from quantarhei import Molecule
from quantarhei import Mode

from quantarhei import CorrelationFunction
from quantarhei import energy_units
from quantarhei import TimeAxis

from quantarhei.qm import ReducedDensityMatrix

import quantarhei as qr

class AggregateTest(unittest.TestCase):
    """Tests for the Manager class
    
    
    """
    
    def setUp(self):
        
        m1 = Molecule(name="Molecule 1", elenergies=[0.0, 1.0])
        m2 = Molecule(name="Molecule 2", elenergies=[0.0, 1.0])
        
        time = TimeAxis(0.0, 1000, 1.0)
        params = dict(ftype="OverdampedBrownian", reorg=20, cortime=100, T=300)
        with energy_units("1/cm"):
            fc = CorrelationFunction(time, params)
            
        m1.set_transition_environment((0,1), fc)
        m1.position = [0.0, 0.0, 0.0]
        m1.set_dipole(0,1,[10.0, 0.0, 0.0])
        m2.set_transition_environment((0,1), fc)
        m2.position = [10.0, 0.0, 0.0]
        m2.set_dipole(0,1,[10.0, 0.0, 0.0])
        
        self.agg = Aggregate(name="TestAgg", molecules=[m1,m2])
        
        self.agg.set_coupling_by_dipole_dipole()
        
        self.agg.build()
        
        
        m3 = Molecule(name="Molecule 1", elenergies=[0.0, 1.0])
        m4 = Molecule(name="Molecule 2", elenergies=[0.0, 1.0])     
        m3.add_Mode(Mode(0.01))
        m4.add_Mode(Mode(0.01))
        
        mod3 = m3.get_Mode(0)
        mod4 = m4.get_Mode(0)
        
        mod3.set_nmax(0, 4)
        mod3.set_nmax(1, 4)
        mod3.set_HR(1, 0.1)
        mod4.set_nmax(0, 4)
        mod4.set_nmax(1, 4)
        mod4.set_HR(1, 0.3)
        
        self.vagg = Aggregate(molecules=[m3, m4])
        self.vagg.build()


#########################################################################
#
#                               TESTS
#
#########################################################################

    def test_saving_of_aggregate(self):
        """(Aggregate) Testing its saving capability 
        
        
        """
        use_temporary_file = True
        
        if use_temporary_file: 
            
            #drv = "core"
            #bcs = False
        
            #with h5py.File('tempfile.hdf5', 
            #               driver=drv, 
            #               backing_store=bcs) as f:
            with tempfile.TemporaryFile() as f:
                                             
                self.agg.save(f, test=True) #,report_unsaved=True)
                
                # reread it
                agg = Aggregate()
                agg = agg.load(f, test=True)

        else:

            #with h5py.File('tempfile.hdf5') as f:
            with open('tempfile.qrp','wb') as f:                                           
                self.agg.save(f)
            
            #with h5py.File('tempfile.hdf5') as f:
            with open('tempfile.qrp','rb') as f:
                agg = Aggregate()
                agg.load(f)
            
            
        
        self.assertEqual(self.agg.nmono, agg.nmono)
        self.assertEqual(self.agg.mult, agg.mult)
        self.assertEqual(self.agg.sbi_mult, agg.sbi_mult)
        self.assertEqual(self.agg.name, agg.name)
        self.assertEqual(self.agg._has_egcf_matrix, agg._has_egcf_matrix)
        self.assertEqual(self.agg._has_system_bath_interaction,
                         agg._has_system_bath_interaction)
        self.assertEqual(self.agg.coupling_initiated, agg.coupling_initiated)
        self.assertEqual(self.agg._has_relaxation_tensor,
                         agg._has_relaxation_tensor)
        self.assertEqual(self.agg._relaxation_theory, agg._relaxation_theory)
        self.assertEqual(self.agg._built, agg._built)
        self.assertEqual(self.agg.Nel,agg.Nel)
        self.assertEqual(self.agg.Ntot,agg.Ntot)
        numpy.testing.assert_array_equal(self.agg.Nb,agg.Nb)
        for key in agg.mnames.keys():
            self.assertEqual(self.agg.mnames[key],agg.mnames[key])
        numpy.testing.assert_array_equal(self.agg.resonance_coupling, 
                                         agg.resonance_coupling) 
        mnames = ["Molecule 1", "Molecule 2"]
        for k in range(agg.nmono):
            m = agg.monomers[k]
            self.assertIn(m.name, mnames)
            mnames.remove(m.name)
            

    def test_trace_over_vibrations(self):
        """(Aggregate) Testing trace over vibrational DOF
        
        """
        
        agg = self.vagg
        
        ham = agg.get_Hamiltonian()
        
        rho = ReducedDensityMatrix(dim=ham.dim)
        rho._data[5, 5] = 1.0
        redr = agg.trace_over_vibrations(rho)
        
        numpy.testing.assert_almost_equal(numpy.trace(redr._data), 1.0,
                                          decimal=7)
        
        N = agg.Nb[0]
        
        rho = ReducedDensityMatrix(dim=ham.dim)
        rho._data[N+5, N+5] = 1.0
        redr = agg.trace_over_vibrations(rho)       

        numpy.testing.assert_almost_equal(numpy.trace(redr._data), 1.0,
                                          decimal=3)


        N = agg.Nb[1]
        
        rho = ReducedDensityMatrix(dim=ham.dim)
        rho._data[N+5, N+5] = 1.0
        redr = agg.trace_over_vibrations(rho)       

        numpy.testing.assert_almost_equal(numpy.trace(redr._data), 1.0,
                                          decimal=2)

        rho = ReducedDensityMatrix(dim=ham.dim)
        n_part = 1./(agg.Nb[0])
        for state_ind in range(agg.Nb[0]):
            rho._data[state_ind, state_ind] = n_part

        redr = agg.trace_over_vibrations(rho)  
        numpy.testing.assert_almost_equal(numpy.trace(redr._data), 1.0,
                                    decimal=7)

        
    def test_get_Density_Matrix_thermal(self):
        """(Aggregate) Testing the get_Densitymatrix method with `thermal` condition type
        
        """
        from quantarhei.core.units import kB_intK
        
        vagg = self.vagg
        H = vagg.get_Hamiltonian()
        N = H.data.shape[0]
        
        # zero temperature
        rho0 = vagg.get_DensityMatrix(condition_type="thermal")
        tst0 = numpy.zeros((N,N),dtype=qr.COMPLEX)
        tst0[0,0] = 1.0
        
        numpy.testing.assert_almost_equal(rho0.data, tst0)
        
        # room temperature
        rho0 = vagg.get_DensityMatrix(condition_type="thermal", 
                                      temperature=300)
        # cannonical balance between neighboring states is tested
        kbT = kB_intK*300.0
        for k in range(N-1):
            E2 = H.data[k+1,k+1]
            E1 = H.data[k,k]
            ratio = numpy.exp(-(E2-E1)/kbT)
            numpy.testing.assert_almost_equal(ratio, 
                                       rho0.data[k+1,k+1]/rho0.data[k,k])

            
    def test_get_Density_Matrix_thermal_excited(self):
        """(Aggregate) Testing the get_Densitymatrix method with `thermal_excited_state` condition type
        
        """
        from quantarhei.core.units import kB_intK
        
        vagg = self.vagg
        H = vagg.get_Hamiltonian()
        N = H.data.shape[0]
        Nb0 = vagg.Nb[0]
           
        # zero temperature
        rho0 = vagg.get_DensityMatrix(condition_type="thermal_excited_state")
        tst0 = numpy.zeros((N,N),dtype=qr.COMPLEX)
        tst0[Nb0,Nb0] = 1.0        
        
        numpy.testing.assert_almost_equal(rho0.data, tst0)
        
        # room temperature
                                      
        kbT = kB_intK*300.0
        
        H = self.agg.get_Hamiltonian()
        N = H.data.shape[0]
        Nb0 = self.agg.Nb[0]
        
        with qr.eigenbasis_of(H):
            rho0 = self.agg.get_DensityMatrix(condition_type="thermal_excited_state", 
                                              temperature=300)
       
        tst0 = numpy.zeros((N,N),dtype=qr.COMPLEX)
        
        with qr.eigenbasis_of(H):
            ssum = 0.0
            for i in range(Nb0, N):
                tst0[i,i] = numpy.exp(-(H.data[i,i]-H.data[Nb0,Nb0])/kbT)
                ssum += tst0[i,i]
            tst0 = tst0/ssum        
            numpy.testing.assert_almost_equal(rho0.data, tst0)
        
        
        
#        ssum = 0.0
#        for i in range(N):
#            tst0[i,i] = numpy.exp(-H.data[i,i]/kbT)
#            ssum += tst0[i,i]
#        tst0 = tst0/ssum
#        
#        DD = vagg.TrDMOp.data
#        # abs value of the transition dipole moment
#        dabs = numpy.sqrt(DD[:,:,0]**2 + \
#                          DD[:,:,1]**2 + DD[:,:,2]**2)
#        # excitation from bra and ket
#        tst0 = numpy.dot(dabs, numpy.dot(tst0,dabs)) 
#        
#        numpy.testing.assert_almost_equal(rho0.data, tst0)
        
        
        

    def test_get_temperature(self):
        """(Aggregate) Testing temperature retrieval 
        
        
        """
        agg = self.agg
        
        T = agg.get_temperature()
        
        numpy.testing.assert_almost_equal(T, 300.0)
        
  
    def test_init_coupling_matrix(self):
        """(Aggregate) Testing coupling matrix initialization 
        
        """
        agg = TestAggregate(name="dimer-2-env")
        
        self.assertFalse(agg.coupling_initiated)
        
        agg.init_coupling_matrix()
        
        self.assertTrue(agg.coupling_initiated)
        self.assertAlmostEqual(0.0, agg.resonance_coupling[0,1])

     
    def test_set_resonance_coupling(self):        
        """(Aggregate) Testing resonance coupling setting with different units
        
        """
        agg = TestAggregate(name="dimer-2-env")
        agg.init_coupling_matrix()        
        
        agg.set_resonance_coupling(0, 1, 1.0)
        self.assertAlmostEqual(1.0, agg.resonance_coupling[0,1])
        self.assertAlmostEqual(1.0, agg.resonance_coupling[1,0])
        
        with qr.energy_units("1/cm"):
            agg.set_resonance_coupling(0, 1, 100.0)
        self.assertAlmostEqual(qr.convert(100.0, "1/cm", "int"),
                               agg.resonance_coupling[0,1])
        self.assertAlmostEqual(qr.convert(100.0, "1/cm", "int"),
                               agg.resonance_coupling[1,0])
        
     
    def test_get_resonance_coupling(self):
        """(Aggregate) Testing resonance coupling retrieval in different units
        
        """
        agg = TestAggregate(name="dimer-2-env")
        agg.init_coupling_matrix()        
        agg.set_resonance_coupling(0, 1, 1.0)

        coup = agg.get_resonance_coupling(1,0)
        self.assertAlmostEqual(coup, 1.0)

        with qr.energy_units("1/cm"):
            coup = agg.get_resonance_coupling(1,0)
        self.assertAlmostEqual(coup, qr.convert(1.0, "int", "1/cm"))
        
        with qr.energy_units("eV"):
            coup = agg.get_resonance_coupling(1,0)
        self.assertAlmostEqual(coup, qr.convert(1.0, "int", "eV"))
        
        with qr.energy_units("THz"):
            coup = agg.get_resonance_coupling(1,0)
        self.assertAlmostEqual(coup, qr.convert(1.0, "int", "THz"))
        
       
    def test_set_resonance_coupling_matrix(self):
        """(Aggregate) Testing setting the whole resonance coupling matrix
        
        """
        agg = TestAggregate(name="dimer-2-env")

        mat = [[0.0, 1.0],
               [1.0, 0.0]]
        
        agg.set_resonance_coupling_matrix(mat)
        self.assertAlmostEqual(agg.resonance_coupling[1,0], 1.0)
        self.assertAlmostEqual(agg.resonance_coupling[0,1], 1.0)
        self.assertAlmostEqual(agg.resonance_coupling[0,0], 0.0)
        self.assertAlmostEqual(agg.resonance_coupling[1,1], 0.0)

        mat = ((0.0, 500.0),
               (500.0, 0.0))
        
        with qr.energy_units("1/cm"):
            agg.set_resonance_coupling_matrix(mat)

        inint = qr.convert(500, "1/cm", "int")
        self.assertAlmostEqual(agg.resonance_coupling[1,0], inint)
        self.assertAlmostEqual(agg.resonance_coupling[0,1], inint)
        self.assertAlmostEqual(agg.resonance_coupling[0,0], 0.0)
        self.assertAlmostEqual(agg.resonance_coupling[1,1], 0.0)

    
    def test_dipole_dipole_coupling(self):
        """(Aggregate) Testing dipole-dipole coupling calculation

        """
        agg = TestAggregate(name="dimer-2-env")

        coup1 = agg.dipole_dipole_coupling(0,1)
        coup2 = agg.dipole_dipole_coupling(0,1, epsr=2.0)
        self.assertAlmostEqual(coup1, coup2*2.0)
        
        with qr.energy_units("1/cm"):
            coup = agg.dipole_dipole_coupling(0,1)
        self.assertAlmostEqual(qr.convert(coup1,"int","1/cm"), coup)


    def test_set_coupling_by_dipole_dipole(self):
        """(Aggregate) Testing dipole-dipole coupling calculation for the whole agregate
        
        """
        agg = TestAggregate(name="dimer-2-env")
        agg.set_coupling_by_dipole_dipole()

        coup1 = agg.dipole_dipole_coupling(0,1)
        self.assertAlmostEqual(coup1, agg.resonance_coupling[0,1])
        
        with qr.energy_units("1/cm"):
            agg.set_coupling_by_dipole_dipole()
            
        coup1 = agg.dipole_dipole_coupling(0,1)
        self.assertAlmostEqual(coup1, agg.resonance_coupling[0,1])


    def test_calculate_resonance_coupling(self):
        """(Aggregate) Testing resonance coupling calculation by various methods
        
        """
        agg = TestAggregate(name="dimer-2-env")
        
        agg.calculate_resonance_coupling(method="dipole-dipole")        
        coup1 = agg.dipole_dipole_coupling(0,1)
        
        self.assertEqual(coup1, agg.resonance_coupling[1,0])
        
        with qr.energy_units("1/cm"):
            agg.calculate_resonance_coupling(method="dipole-dipole")        

        self.assertEqual(coup1, agg.resonance_coupling[1,0])
            
            
    def test_add_Molecule(self):
        """(Aggregate) Testing add_Molecule() method
        
        """
        agg = TestAggregate(name="dimer-2-env")   
        
        mol = qr.Molecule(elenergies=[0.0, 1.0])
        
        nmols1 = agg.nmono
        agg.add_Molecule(mol)
        nmols2 = agg.nmono
        
        self.assertEqual(nmols1 + 1, nmols2)
        self.assertEqual(len(agg.monomers), nmols2)
        self.assertLessEqual(len(agg.mnames), nmols2)
        
        
        
    def test_projection_operator_generation(self):
        """(Aggregate) Testing projector operator generation
        
        """        
        agg = self.vagg
        
        T1 = 800.0
        T2 = 1500.0
        
        # figure out which sites have vibrations
        
        # what is the position of my site (verify that site has vibrations)
        
        # generate all vibrational signatures left and right
        orates = [1.0/T1, 1.0/T2]
        sbi = qr.qm.SystemBathInteraction(osites=[0,1], orates=orates,
                                          system=agg)
        # 
        from quantarhei.qm.liouvillespace.lindbladform import \
        VibrationalDecayLindbladForm
        
        ham = agg.get_Hamiltonian()
        
        vlf = VibrationalDecayLindbladForm(ham, sbi)
        
        self.assertEqual(T1, numpy.sqrt(3.0)/(2.0*numpy.max(vlf.Lm[2,:,:])))
        self.assertEqual(T2, numpy.sqrt(3.0)/(2.0*numpy.max(vlf.Lm[5,:,:])))        
        
        