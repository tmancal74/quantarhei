# -*- coding: utf-8 -*-
"""
    Representation of Hierarchical Equations of Motion


"""
import numpy
from ... import REAL, COMPLEX
from ..hilbertspace.operators import ReducedDensityMatrix

class KTHierarchy:
    """ Kubo-Tanimura Hierarchy
    
    
    Parameters
    ----------
    
    ham : Hamiltonian
        System Hamiltonian
        
    sbi : SystemBathInteraction
        System bath interaction object
        
    depth : int
        Depth of the hierarchy
        
        
    >>> import quantarhei as qr
    >>> m1 = qr.Molecule([0.0, 1.0])
    >>> m2 = qr.Molecule([0.0, 1.0])
    >>> agg = qr.Aggregate([m1, m2])
    >>> agg.set_resonance_coupling(0,1,0.1)
    >>> agg.build()
    >>> ham = agg.get_Hamiltonian()
    >>> print(ham)
    <BLANKLINE>
    quantarhei.Hamiltonian object
    =============================
    units of energy 1/fs
    Rotating Wave Approximation (RWA) enabled : True
    Number of blocks : 2
    Block average energies:
     0 : 0.0
     1 : 1.0
    data = 
    [[ 0.   0.   0. ]
     [ 0.   1.   0.1]
     [ 0.   0.1  1. ]]
    >>> sbi = qr.qm.TestSystemBathInteraction("dimer-2-env")
    >>> Hy = KTHierarchy(ham, sbi, 4)
    >>> print(Hy.dim)
    3
    >>> print(Hy.nbath)
    2
    >>> print(Hy.gamma)
    [ 0.01  0.01]
    >>> print(Hy.lam)
    [ 0.0037673  0.0037673]
    >>> print(Hy.depth)
    4
    >>> print(Hy.temp)
    300
    >>> print(Hy.hsize)
    15
    >>> print(Hy.indxs)
    [[[0, 0]], [[1, 0], [0, 1]], [[2, 0], [1, 1], [0, 2]], [[3, 0], [2, 1], [1, 2], [0, 3]], [[4, 0], [3, 1], [2, 2], [1, 3], [0, 4]]]
    >>> print(Hy.Vs[0,:,:])
    [[ 0.  0.  0.]
     [ 0.  1.  0.]
     [ 0.  0.  0.]]
    >>> print(Hy.rho)
    <BLANKLINE>
    quantarhei.ReducedDensityMatrix object
    ======================================
    data = 
    [[ 0.+0.j  0.+0.j  0.+0.j]
     [ 0.+0.j  0.+0.j  0.+0.j]
     [ 0.+0.j  0.+0.j  0.+0.j]]


    """
    
    def __init__(self, ham, sbi, depth=2):
        
        self.ham = ham
        self.sbi = sbi
        self.depth = depth
        
        
        # dimension of the ADOs
        self.dim = ham.dim
        
        # number of baths
        self.nbath = self.sbi.N
        
        self.gamma = numpy.zeros(self.nbath, dtype=REAL)
        for ii in range(self.nbath):
            self.gamma[ii] = 1.0/self.sbi.get_correlation_time(ii)
            
        self.lam = numpy.zeros(self.nbath, dtype=REAL)
        for ii in range(self.nbath):
            self.lam[ii] = self.sbi.get_reorganization_energy(ii)
            
        self.temp = self.sbi.get_temperature()
        
        # generation of hierarchy indices
        self.indxs = self._generate_indices(self.nbath, level=self.depth)
        
        self.hsize = 0
        for levels in self.indxs:    
            self.hsize += len(levels)
        
        self.Vs = self.sbi.KK
        
        self.rho = ReducedDensityMatrix(data=numpy.zeros((self.dim, self.dim),
                                                         dtype=COMPLEX))
        # This needs to be basis controlled
        self.ado = None
        self.reset_ados()


    def _generate_indices(self, N, level=0):
        """Generation of hierarchy indices for a give problem
        
        """
        if N != 2:
            raise Exception("Experimental code, N different from 2"+
                            " not implemented")
        if level > 4:
            raise Exception("Experimental code, level > 4 not implemented")
            
        lret = []
        
        level0 = []
        level0.append([0,0])
        lret.append(level0)
        if level == 0:
            return lret
        level1 = []
        level1.append([1,0])
        level1.append([0,1])
        lret.append(level1)
        if level == 1:
            return lret
        level2 = []
        level2.append([2,0])
        level2.append([1,1])
        level2.append([0,2])
        lret.append(level2)
        if level == 2:
            return lret
        level3 = []
        level3.append([3,0])
        level3.append([2,1])
        level3.append([1,2])
        level3.append([0,3])
        lret.append(level3)
        if level == 3:
            return lret
        level4 = []
        level4.append([4,0])
        level4.append([3,1])
        level4.append([2,2])
        level4.append([1,3])        
        level4.append([0,4])
        lret.append(level4)                
        return lret
        
    
    def set_rdo(self, rdo):
        """Sets the density operator of the hierarchy
        """
        self.rho = rdo
        
    def reset_ados(self):
        """Creates memory of ADOs and sets them to zero
        
        """
        self.ado = numpy.zeros((self.hsize, self.dim, self.dim-1),
                               dtype=COMPLEX)
        
        
class KTHierarchyPropagator:
    
    def __init__(self, timeaxis, hierarchy):
        
        pass
    
    
    def propagate(self, rhoi):
        
        rho = rhoi
        
        return rho
    