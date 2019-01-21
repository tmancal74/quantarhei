# -*- coding: utf-8 -*-
"""
    Representation of Hierarchical Equations of Motion


"""
import numpy
from ... import REAL, COMPLEX
from ..hilbertspace.operators import ReducedDensityMatrix
from ..propagators.dmevolution import DensityMatrixEvolution

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
    >>> print(Hy.levels)
    [ 0  1  3  6 10]
    >>> print(Hy.levlengths)
    [1 2 3 4 5]
    >>> print(Hy.hinds)
    [[0 0]
     [1 0]
     [0 1]
     [2 0]
     [1 1]
     [0 2]
     [3 0]
     [2 1]
     [1 2]
     [0 3]
     [4 0]
     [3 1]
     [2 2]
     [1 3]
     [0 4]]

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
        self.kBT = 1.0
        
        # generation of hierarchy indices
        indxs = self._generate_indices(self.nbath, level=self.depth)
        
        self.hsize = 0
        for levels in indxs:    
            self.hsize += len(levels)
        
        self.Vs = self.sbi.KK
        
        self.rho = ReducedDensityMatrix(data=numpy.zeros((self.dim, self.dim),
                                                         dtype=COMPLEX))
        # This needs to be basis controlled
        self.ado = None
        self.reset_ados()
        
        #
        # numpy representation of the hierarchy indices
        #
        
        # indices where the levels start
        self.levels = numpy.zeros(depth+1, dtype=numpy.int)
        # lengths of the levels
        self.levlengths = numpy.zeros(depth+1, dtype=numpy.int)
        # indices
        self.hinds = self._convert_2_matrix(indxs)
        
        self.nm1 = numpy.zeros((self.hsize, self.nbath), dtype=numpy.int)
        self.np1 = numpy.zeros((self.hsize, self.nbath), dtype=numpy.int)
        self._make_nmp1()
        
        self.Gamma = numpy.zeros(self.hsize, dtype=REAL)
        self._make_Gamma()
        


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
      
    def _convert_2_matrix(self, indxs):
        hsize = 0
        lvl = 0
        start = 0
        self.levels[lvl] = start
        for levels in indxs:
            if lvl > 0:
                self.levels[lvl] = start
            lngth = len(levels)
            self.levlengths[lvl] = lngth
            hsize += lngth
            lvl +=1
            start = start+lngth
            
        mat = numpy.zeros((hsize, self.nbath), dtype=numpy.int)
        ii = 0
        for level in indxs:
            for inds in level:
                kk = 0
                for ind in inds:
                    mat[ii, kk] = ind
                    kk += 1
                ii += 1
        
        return mat
        
    
    def _make_nmp1(self):
        """ Makes the list of indices obtained from n by -1 or +1 operations
        
        """
        
        for nn in range(self.hsize):
            for kk in range(self.nbath):
                indxm = numpy.zeros(self.nbath, dtype=numpy.int)
                indxm[:] = self.hinds[nn,:]
                indxm[kk] -= 1

                indxp = numpy.zeros(self.nbath, dtype=numpy.int)
                indxp[:] = self.hinds[nn,:]
                indxp[kk] += 1
                
                venm = -1
                venp = -1
                for ll in range(nn):
                    if numpy.array_equal(self.hinds[ll,:], indxm):
                        venm = ll
                    if numpy.array_equal(self.hinds[ll,:], indxp):
                        venp = ll
                
                self.nm1[nn, kk] = venm
                self.np1[nn, kk] = venp
        
        
    def _make_Gamma(self):
        """ Decay factor of a given ADO
        
        """
        
        for nn in range(self.hsize):
            for kk in range(self.nbath):
                self.Gamma[nn] += self.hinds[nn,kk]*self.gamma[kk]

        
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
    """Propagator of the Kubo-Tanimura hierarchy
    
    >>> import quantarhei as qr
    >>> m1 = qr.Molecule([0.0, 1.0])
    >>> m2 = qr.Molecule([0.0, 1.0])
    >>> agg = qr.Aggregate([m1, m2])
    >>> agg.set_resonance_coupling(0,1,0.1)
    >>> agg.build()
    >>> ham = agg.get_Hamiltonian()
    >>> sbi = qr.qm.TestSystemBathInteraction("dimer-2-env")
    >>> Hy = KTHierarchy(ham, sbi, 4)
    >>> rhoi = qr.ReducedDensityMatrix(dim=ham.dim)
    >>> rhoi.data[2,2] = 1.0
    >>> print(rhoi)
    <BLANKLINE>
    quantarhei.ReducedDensityMatrix object
    ======================================
    data = 
    [[ 0.+0.j  0.+0.j  0.+0.j]
     [ 0.+0.j  0.+0.j  0.+0.j]
     [ 0.+0.j  0.+0.j  1.+0.j]]
    >>> time = qr.TimeAxis(0.0, 1000, 1.0)
    >>> kprop = KTHierarchyPropagator(time, Hy)
    >>> rhot = kprop.propagate(rhoi)
    >>> #print(rhot.data[:,2,2])
    
    
    """
    
    def __init__(self, timeaxis, hierarchy):
        
        self.timeaxis = timeaxis
        self.Nt = timeaxis.length
        self.dt = timeaxis.step
        self.hy = hierarchy
        
    
    
    def propagate(self, rhoi, L=4):
        """Propagates the Kubo-Tanimura Hierarchy including the RDO
        
        """
        rhot = DensityMatrixEvolution(timeaxis=self.timeaxis, rhoi=rhoi)

        ado1 = self.hy.ado
        ado2 = self.hy.ado

        self.Nref = 1

        indx = 1
        for ii in self.timeaxis.data[1:self.Nt]:

            for jj in range(0,self.Nref):

                for ll in range(1,L+1):

                    ado1 = self._ado_cros_rhs(ado1, (self.dt/ll)) \
                         + self._ado_self_rhs(ado1, (self.dt/ll))

                    ado2 = ado2 + ado1      
                ado1 = ado2

            self.hy.ado = ado2
            rhot.data[indx,:,:] = ado2[0,:,:]                        
            indx += 1             

        return rhot


    def _ado_self_rhs(self, ado1, dt):
        """Self contribution of the equation for the hierarchy ADOs

        """
        ado3 = numpy.zeros(ado1.shape, dtype=ado1.dtype)
        HH = self.hy.ham.data
        
        for nn in range(self.hy.hsize):
            
            # 
            ado3[nn,:,:] = -dt*(1j*(numpy.dot(HH, ado1[nn,:,:])
                                  - numpy.dot(ado1[nn,:,:],HH)) 
                                + self.hy.Gamma[nn]*ado1[nn,:,:])
           
                           
        return ado3

    
    
    def _ado_cros_rhs(self, ado1, dt):
        """All cross-terms of the Hierarchy 
        
        """
        
        ado3 = numpy.zeros(ado1.shape, dtype=ado1.dtype)
        rl = numpy.zeros((ado1.shape[1],ado1.shape[2]), dtype=ado1.dtype)
        rr = numpy.zeros((ado1.shape[1],ado1.shape[2]), dtype=ado1.dtype)
        
        
        for nn in range(self.hy.hsize):
            for kk in range(self.hy.nbath):
                
                # Theta+ and Psi+
                nk = self.hy.hinds[nn,kk]   
                jj = self.hy.nm1[nn,kk]
                
                if nk*jj > 0:
                    
                    rr = numpy.dot(self.hy.Vs[kk,:,:], ado1[jj,:,:])
                    rl = numpy.dot(ado1[jj,:,:],self.hy.Vs[kk,:,:])
                    
                    # Theta
                    ado3[nn,:,:] += dt*nk*self.hy.lam[kk]*self.hy.gamma[kk]* \
                                    (rr-rl)
                        
                            
                    # Psi
                    ado3[nn,:,:] += (1j*dt)*2.0*nk*self.hy.lam[kk]*self.kBT* \
                                    (rr+rl)
                    
                # Psi-
                jj = self.hy.np1[nn,kk]
                
                if jj >= 0:

                    rr = numpy.dot(self.hy.Vs[kk,:,:], ado1[jj, :,:])
                    rl = numpy.dot(ado1[jj,:,:],self.hy.Vs[kk,:,:])
                    
                    ado3[nn,:,:] += (1j*dt)*(rr-rl)
                    
        return ado3