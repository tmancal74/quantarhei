# -*- coding: utf-8 -*-
"""
    Representation of Hierarchical Equations of Motion


"""
import numpy
from ... import REAL, COMPLEX
from ..propagators.dmevolution import DensityMatrixEvolution
from ..hilbertspace.operators import ReducedDensityMatrix
from ..hilbertspace.operators import UnityOperator
from ...core.units import kB_int
from ..corfunctions.correlationfunctions import CorrelationFunction

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
    Warning: OverdampedBrownian spectral density
     - only high-temperature limit of this function is used in HEOM
    Warning: OverdampedBrownian spectral density
     - only high-temperature limit of this function is used in HEOM
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

        # FIXME
        # check that sbi only has correlation functions of Lorentz type
        for ii in range(self.nbath):
            cc = self.sbi.CC.get_correlation_function(ii,ii)
            prms = cc.params[0]
            if prms["ftype"] == CorrelationFunction.allowed_types[1]:
                tp = CorrelationFunction.allowed_types[1]
                print("Warning: "+tp+" spectral density\n"+
                      " - only high-temperature limit "+
                      "of this function is used in HEOM")
            elif prms["ftype"] in [CorrelationFunction.allowed_types[0],
                                   "Lorentz-Drude"]:
                pass
            else:    
                raise Exception("Spectral density/Correlation function type:"+
                                prms["ftype"]+
                                "\nHEOM is not implemented for this function")
        
        self.gamma = numpy.zeros(self.nbath, dtype=REAL)
        for ii in range(self.nbath):
            self.gamma[ii] = 1.0/self.sbi.get_correlation_time(ii)
            
        self.lam = numpy.zeros(self.nbath, dtype=REAL)
        for ii in range(self.nbath):
            self.lam[ii] = self.sbi.get_reorganization_energy(ii)
            
        self.temp = self.sbi.get_temperature()
        self.kBT = self.temp*kB_int
        
        # generation of hierarchy indices
        indxs = self.generate_indices(self.nbath, level=self.depth)
        
        self.hsize = 0
        for levels in indxs:    
            self.hsize += len(levels)
        
        self.Vs = self.sbi.KK
        
        #self.rho = ReducedDensityMatrix(data=numpy.zeros((self.dim, self.dim),
        #                                                 dtype=COMPLEX))
        
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
        
        self.hpop = None
        

    
    def generate_indices(self, N, level=0):
        """Generates indices of the hierarchy up a certain level
        
        
        Parameters
        ----------
        
        N : int
            Number of indices in the hierarchy
            
        level : int
            Highest level of the hierarchy to be generated
        
        """
        
        if False:
            return self._generate_indices_2_0to4(N, level)
        else:
            
            lret = []
            
            if False:
                level_prev = []
                inilist = [0]*N          # lowest level 
                level_prev.append(inilist)
                lret.append(level_prev)            
                
                for kk in range(level):
                    last_level = kk
                    new_level_prev = []
                    for old_level in level_prev:
                        doit = False
                        for nn in range(N):
                            if old_level[nn] == last_level:
                                doit = True
                            if doit:
                                nlist = old_level.copy()
                                nlist[nn] += 1
                                new_level_prev.append(nlist)
                
                    level_prev = new_level_prev
                    lret.append(level_prev)
            else:
                
                level_prev = []
                inilist = [0]*N
                level_prev.append(inilist)
                lret.append(level_prev)
                
                for kk in range(level):
                    last_level = kk
                    new_level_prev = []
                    for old_level in level_prev:
                        for nn in range(N):
                            nlist = old_level.copy()
                            nlist[nn] += 1
                            #check if it is already in
                            if nlist not in new_level_prev:
                                new_level_prev.append(nlist)
                                
                    level_prev = new_level_prev
                    lret.append(level_prev)
                
            return lret
            
            


    def _generate_indices_2_0to4(self, N=2, level=0):
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
        """Convert the list of levels of the hierarchy into an numpy array
        
        
        """
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
                for ll in range(nn):
                    if numpy.array_equal(self.hinds[ll,:], indxm):
                        venm = ll
                venp = -1
                for ll in range(self.hsize):
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

   
    def reset_ados(self):
        """Creates memory of ADOs and sets them to zero
        
        """
        self.ado = numpy.zeros((self.hsize, self.dim, self.dim),
                               dtype=COMPLEX)


    def get_kernel(self, timeaxis):
        """Returns integration kernel for the time-non-local equation
        
        """
        N = self.dim
        Nt = timeaxis.length
        kernel = numpy.zeros((Nt,N,N,N,N), dtype=COMPLEX)
        
        khprop = KTHierarchyPropagator(timeaxis, self)
        
        rhoi = ReducedDensityMatrix(dim=self.dim)
        for ii in range(N):
            for jj in range(N):
                print("Starting:", ii, jj)
                rhoi = ReducedDensityMatrix(dim=self.dim)
                rhoi.data[ii,jj] = 1.0
                
                rhot = khprop.propagate(rhoi, free_hierarchy=True)
                kernel[:,:,:,ii,jj] = -rhot.data[:,:,:]

                print("... finished.")
                
        qhp = self._QHPsop()
        phq = self._PHQsop()
        for tk in range(Nt):
            k1 = numpy.tensordot(kernel[tk,:,:,:,:], qhp)
            kernel[tk,:,:,:,:] = numpy.tensordot(phq, k1)
        
        return kernel


    def _QHPsop(self):
        
        N = self.dim
        delta = UnityOperator(dim=N).data
        qhp = numpy.zeros((N, N, N, N), dtype=COMPLEX)
        
        for kk in range(self.nbath):
                
            # Theta+ and Psi+
            nk = self.hinds[1,kk]
            jj = self.nm1[1,kk]
            
            print(kk, nk, jj)

            if nk*jj > 0:

                print("Calculating QHP:")
                for ii_i in range(N):
                    for jj_i in range(N):
                        for kk_i in range(N):
                            for ll_i in range(N):
                                # Theta
                                qhp[ii_i,jj_i,kk_i,ll_i] += \
                                nk*self.lam[kk]*self.gamma[kk]*\
                                (self.Vs[kk,ii_i,kk_i]*delta[jj_i,ll_i]
                                 + delta[kk_i,ii_i]*self.Vs[kk,ll_i,jj_i])
                                # Psi
                                qhp[ii_i,jj_i,kk_i,ll_i] += \
                                1j*2.0*nk*self.lam[kk]*self.kBT* \
                                (self.Vs[kk,ii_i,kk_i]*delta[jj_i,ll_i]
                                 - delta[kk_i,ii_i]*self.Vs[kk,ll_i,jj_i])
                         
                print(" ...done")
                
        return qhp                


    def _PHQsop(self):
        
        N = self.dim
        delta = UnityOperator(dim=N).data
        phq = numpy.zeros((N, N, N, N), dtype=COMPLEX)
        for ii in range(N):
            for jj in range(N):
                for kk in range(N):
                    for ll in range(N):
                        phq[ii,jj,kk,ll] = delta[ii,kk]*delta[jj,ll]
                    
        return phq        


        
class KTHierarchyPropagator:
    """Propagator of the Kubo-Tanimura hierarchy
    
    >>> import numpy
    >>> import quantarhei as qr
    >>> with qr.energy_units("1/cm"):
    ...     m1 = qr.Molecule([0.0, 10000.0])
    ...     m2 = qr.Molecule([0.0, 10000.0])
    ...     #m3 = qr.Molecule([0.0, 10000.0])
    >>> agg = qr.Aggregate([m1, m2])
    >>> with qr.energy_units("1/cm"):
    ...     agg.set_resonance_coupling(0,1,80.0)
    ...     #agg.set_resonance_coupling(0,2,100.0)
    >>> agg.build()
    >>> ham = agg.get_Hamiltonian()
    >>> sbi = qr.qm.TestSystemBathInteraction("dimer-2-env")
    >>> Hy = KTHierarchy(ham, sbi, 4)
    Warning: OverdampedBrownian spectral density
     - only high-temperature limit of this function is used in HEOM
    Warning: OverdampedBrownian spectral density
     - only high-temperature limit of this function is used in HEOM
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
    
    #>>> import matplotlib.pyplot as plt
    #>>> N = time.length
    #>>> with qr.eigenbasis_of(ham):
    #...     plt.plot(time.data[0:N], rhot.data[0:N,1,1],"-b")
    #...     plt.plot(time.data[0:N], rhot.data[0:N,2,2],"-r")
    #...     plt.plot(time.data[0:N], numpy.real(rhot.data[0:N,1,2]),"-g")
    #>>> plt.show()

    
    
    """
    
    def __init__(self, timeaxis, hierarchy):
        
        self.timeaxis = timeaxis
        self.Nt = timeaxis.length
        self.dt = timeaxis.step
        self.hy = hierarchy

        #
        # RWA
        #
        if self.hy.ham.has_rwa:
            
            self.RWA = self.hy.ham.rwa_indices
            self.RWU = numpy.zeros(self.RWA.shape, dtype=self.RWA.dtype)
            
            HH = self.hy.ham.data
            shape = HH.shape
            HOmega = numpy.zeros(shape, dtype=REAL)
            for ii in range(shape[0]):
                HOmega[ii,ii] = self.hy.ham.rwa_energies[ii]
                                
            self.HOmega = HOmega        
    
    
    def propagate(self, rhoi, L=4, report_hierarchy=False,
                                   free_hierarchy=False):
        """Propagates the Kubo-Tanimura Hierarchy including the RDO
        
        """
        rhot = DensityMatrixEvolution(timeaxis=self.timeaxis, rhoi=rhoi)
        
        if free_hierarchy:
            
            # first act with lifting superoperators
            self.hy.ado[1,:,:] = rhoi.data
            N = rhoi.dim
            Nt = self.timeaxis.length
            ker = numpy.zeros((Nt,N,N), dtype=COMPLEX)
            ker[0,:,:] = rhoi.data
            
            # Now we propagate normally; slevel is set to 1 so zero's order
            # does not update and stays zero
            slevel = 1
            
        else:
            
            # normally inital condition goes here
            self.hy.ado[0,:,:] = rhoi.data
            slevel = 0
        
        ado1 = self.hy.ado
        ado2 = self.hy.ado

        # no fine time-step for integro-differential solver
        self.Nref = 1
        
        if report_hierarchy:
            # we report population of hierarchy ADO
            Nt = rhot.data.shape[0]
            self.hy.hpop = numpy.zeros((Nt, self.hy.hsize), dtype=REAL)
            for kk in range(self.hy.hsize):
                self.hy.hpop[0,kk] = numpy.trace(self.hy.ado[kk,:,:])

        indx = 1
        for ii in self.timeaxis.data[1:self.Nt]:

            for jj in range(0,self.Nref):

                for ll in range(1,L+1):

                    ado1 = self._ado_cros_rhs(ado1, (self.dt/ll), slevel) \
                         + self._ado_self_rhs(ado1, (self.dt/ll), slevel)

                    ado2 = ado2 + ado1      
                ado1 = ado2

            self.hy.ado = ado2
            
            if free_hierarchy:
                ker[indx,:,:] = ado2[1,:,:] 

            else:
                rhot.data[indx,:,:] = ado2[0,:,:]
                
            if report_hierarchy:
                # we report population of hierarchy ADO
                for kk in range(self.hy.hsize):
                    self.hy.hpop[indx, kk] = numpy.trace(ado2[kk,:,:])                       
            indx += 1             
            
        return rhot


    def _ado_self_rhs(self, ado1, dt, slevel=0):
        """Self contribution of the equation for the hierarchy ADOs

        """
        ado3 = numpy.zeros(ado1.shape, dtype=ado1.dtype)
        
        if self.hy.ham.has_rwa:
            HH = self.hy.ham.data  - self.HOmega
        else:
            HH = self.hy.ham.data
        
        for nn in range(slevel, self.hy.hsize):
            
            # 
            ado3[nn,:,:] = -dt*(1j*(numpy.dot(HH, ado1[nn,:,:])
                                  - numpy.dot(ado1[nn,:,:],HH)) 
                                + self.hy.Gamma[nn]*ado1[nn,:,:])
           
                           
        return ado3

    
    
    def _ado_cros_rhs(self, ado1, dt, slevel=0):
        """All cross-terms of the Hierarchy 
        
        """
        
        ado3 = numpy.zeros(ado1.shape, dtype=ado1.dtype)
        rl = numpy.zeros((ado1.shape[1],ado1.shape[2]), dtype=ado1.dtype)
        rr = numpy.zeros((ado1.shape[1],ado1.shape[2]), dtype=ado1.dtype)
        
        
        for nn in range(slevel, self.hy.hsize):
            for kk in range(self.hy.nbath):
                
                # Theta+ and Psi+
                nk = self.hy.hinds[nn,kk]   
                jj = self.hy.nm1[nn,kk]
                
                
                if nk*jj >= 0:
                    
                    
                    rr = numpy.dot(self.hy.Vs[kk,:,:], ado1[jj,:,:])
                    rl = numpy.dot(ado1[jj,:,:],self.hy.Vs[kk,:,:])
                    
                    # Theta
                    ado3[nn,:,:] += dt*nk*self.hy.lam[kk]*self.hy.gamma[kk]* \
                                    (rr+rl)
                        
                            
                    # Psi
                    ado3[nn,:,:] += (1j*dt)*2.0*nk* \
                                    self.hy.lam[kk]*self.hy.kBT* \
                                    (rr-rl)
                    
                # Psi-
                jj = self.hy.np1[nn,kk]
                
                if jj > 0:
   
                    rr = numpy.dot(self.hy.Vs[kk,:,:], ado1[jj, :,:])
                    rl = numpy.dot(ado1[jj,:,:],self.hy.Vs[kk,:,:])
                    
                    ado3[nn,:,:] += (1j*dt)*(rr-rl)
                    
        return ado3


#    def _QHP(self, rhoi, slevel=1):
#        """One application of the hierarchy operators 
#        
#        """
#        ado1 = numpy.zeros((self.hy.hsize, self.hy.dim, self.hy.dim),
#                           dtype=COMPLEX)
#        ado1[0,:,:] = rhoi.data[:,:]
#        
#        ador =  self._ado_cros_rhs(ado1, dt=1.0, slevel=slevel)
#        ador[0,:,:] = 0.0 # nullify the 0 order so that it cannot contribute
#        
#        return ador



#    def  _PHQ(self, adof, slevel=0):
#        """
#        
#        """
#        return self._ado_cros_rhs(adof, dt=1.0, slevel=slevel)

    
        