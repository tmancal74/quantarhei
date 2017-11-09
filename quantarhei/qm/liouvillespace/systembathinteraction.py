# -*- coding: utf-8 -*-
"""
    Quantarhei package (http://www.github.com/quantarhei)

    systembathinteraction module

"""
import numpy

from ...core.saveable import Saveable
from ...qm.corfunctions.cfmatrix import CorrelationFunctionMatrix

class SystemBathInteraction(Saveable):
    """Describes interaction of an open quantum system with its environment
    
    Stores the system--bath interaction operator in form of a set of operators
    on the Hilbert space of the system and correlation functions of 
    the operator on the bath Hilbert space.
    
    Parameters
    ----------
    
    sys_operators : list
        List of the system part of the system-bath interaction Hamiltonian
        components
    
    bath_correlation_matrix: CorrelationFunctionMatrix 
        Object of the CorrelationFunctionMatrix type holding all correlation
        functions needed for the description of system bath interaction
    
    rates : list/tuple
        List or tuple of rates. The total number of rates has to be 
        the same as the number of operators in the ``sys_operator`` list
        
    system : {Molecule, Aggregate}
        Molecule or Aggregate object in which the system--bath interaction 
        is specified


    """

    def __init__(self, sys_operators=None, bath_correlation_matrix=None,
                 rates=None, system=None):

        # information about aggregate is needed when dealing with 
        # multiple excitons
        self.aggregate = None
        self.molecule = None
        self.system = None
        self.rates = None        
        self.KK = None
        self.CC = None
        self.TimeAxis = None
        self.rates = None
        self.N = 0
        
        # version with bath correlation functions
        if ((sys_operators is not None) 
            and (bath_correlation_matrix is not None)):
            
            # Find the length of the list of operators 
            if isinstance(sys_operators,list):
                self.N = len(sys_operators)
            else:
                raise Exception("sys_operators argument has to a list")
            
            # Second argument has to be a CorrelationFunctionMatrix 
            if not isinstance(bath_correlation_matrix, 
                              CorrelationFunctionMatrix):
                raise Exception("ba_correlation_function argument has to a"+
                                " CorrelationFunctionMatrix")
                
            # Check that sys_operators and bath_correlation matrix has 
            # a compatible number of components
            if bath_correlation_matrix.nob != self.N:
                raise Exception("Incompatile number of bath compoments: " +
                    ("Correlation function matrix - %i vs. operators %i" % 
                    (bath_correlation_matrix.nob,self.N)))
                
            self.TimeAxis = bath_correlation_matrix.timeAxis
            self._set_system(system)

            if self.N > 0:
                
                self._set_operators(sys_operators)
                self.CC = bath_correlation_matrix
                
                
        # version with system-bath operators and rates
        elif ((sys_operators is not None) 
            and (rates is not None)):
            
            # Find the length of the list of operators 
            if isinstance(sys_operators, list):
                self.N = len(sys_operators)
            else:
                raise Exception("First argument has to a list")

            self._set_system(system)
            
            if self.N > 0:
                
                self._set_operators(sys_operators)   
                self.CC = None #bath_correlation_matrix
                
                if len(rates) != self.N:
                    raise Exception("Wrong number of rates specified")
                self.rates = rates


    def _set_system(self, system):
        """Sets the system attribute
        
        """
        from ...builders.aggregates import Aggregate
        from ...builders.molecules import Molecule    
        
        if system is not None:
            if isinstance(system, Aggregate):
                self.aggregate = system
                self.molecule = None
                self.system = self.aggregate

            if isinstance(system, Molecule):
                self.aggregate = None
                self.molecule = system
                self.system = self.molecule


    def _set_operators(self, sys_operators):
        """Sets the system part of the interaction
        
        """
        # First of the operators 
        KK = sys_operators[0]
    
        # Get its dimension 
        dim = KK.data.shape[0]
            
        self.KK = numpy.zeros((self.N, dim, dim), dtype=numpy.float64)
        self.KK[0,:,:] = KK.data       
        
        # Save other operators and check their dimensions 
        for ii in range(1,self.N):
            
            KK = sys_operators[ii]
            
            if True: #isinstance(KK,Operator):
                if dim == KK.data.shape[0]:
                    self.KK[ii,:,:] = KK.data
                else:
                    raise Exception("Operators in the list are" 
                    + " not of the same dimension")
            else:
                raise Exception("sys_operators tuple (the first argument)"
                + " has to contain cu.oqs.hilbertspace.Operator")

        
    def get_coft(self, n, m):
        """Returns bath correlation function corresponding to sites n and m

        
        """
        if self.system is None:
            
            return self.CC.get_coft(n,m)
            
        else:
            
            #FIXME: Molecule needs this method
            bn = self.system.which_band[n]
            bm = self.system.which_band[m]
            
            if ((bn == 0) and (bm == 0)):
                
                #print(bn,"::",n,m)
                return self.CC._cofts[0,:]
                
            elif ((bn == 1) and (bm == 1)):
                #print(bn,"::",n-1,m-1)
                
                return self.CC.get_coft(n-1,m-1)
                
            else:
                
                return self.CC._cofts[0,:]
                
                
                
    def get_coft_elsig(self, n_sig, m_sig):
        """Returns bath correlation based on electronic signatures 


        """ 
                 
        nb = numpy.sum(n_sig)
        mb = numpy.sum(m_sig)
        
        indices = []
        if mb == nb:
            ni = 0
            for na in n_sig:
                mi = 0
                for ma in m_sig:
                    if ((na == 1) and (ma == 1)):
                        indices.append([ni,mi]) 
                    mi += 1
                ni += 1
            
            ret = numpy.zeros((self.TimeAxis.length),dtype=numpy.complex128)
            for ind in indices:
                #print(nb,":",ind[0],ind[1])
                ret += self.get_coft(ind[0],ind[1]) 
                    
                
            return ret            
            
        else:
            return self.CC._cofts[0,:]
