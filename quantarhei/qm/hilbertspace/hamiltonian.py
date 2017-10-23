# -*- coding: utf-8 -*-

from .operators import SelfAdjointOperator
from ...core.managers import BasisManaged
from ...core.managers import EnergyUnitsManaged
from ...utils.types import ManagedRealArray
from .operators import Operator

import numpy


class Hamiltonian(SelfAdjointOperator, BasisManaged, EnergyUnitsManaged):
    """Hamiltonian operator
    
    
    
    """
    
    _has_remainder_coupling = False
    
    data = ManagedRealArray("data")
 
    def __init__(self, dim=None, data=None):
        #self.data = data
#FIXME: how to avoid the Operator breaking the EnergyUnits management ????
        if not ((dim is None) and (data is None)):
            Operator.__init__(self, dim=dim, data=data)
            if not self.check_selfadjoint():
                raise Exception("The data of this operator have"+
                                "to be represented by a selfadjoint matrix")
            
            
    def diagonalize(self, coupling_cutoff=None):
        """Diagonalizes the Hamiltonian matrix 
        
        Parameters
        ----------
        
        coupling_cutoff : float, optional
            Specifies the smallest (absolute) value of coupling 
            which is taken into account. Smaller couplings are removed
            and a remainder coupling matrix is returned together with
            the diagonalization matrix (see Returns section).
            
        Returns
        -------
        
        SS : numpy array
            The diagonalization matrix of the Hamiltonian
            
        JR : numpy array
            Returned only if `coupling_cutoff` is specified. It contains
            couplings that were removed because they are smaller than 
            the cut-off value.
                        
        """
        if coupling_cutoff is None:
            SS = super().diagonalize()
            if self._has_remainder_coupling:
                self.JR = numpy.dot(SS.T,numpy.dot(self.JR,SS))
            self.SS = SS
            return SS
        else:
            self.remove_cutoff_coupling(coupling_cutoff)
            # diagonalize the strong coupling part
            dd,SS = numpy.linalg.eigh(self.data)
            self.data = numpy.zeros(self.data.shape,dtype=numpy.float64)
            for ii in range(0,self.data.shape[0]):
                self.data[ii,ii] = dd[ii]
            # transform the remainder of couling correspondingly
            self.JR = numpy.dot(SS.T,numpy.dot(self.JR,SS))
            self.SS = SS
            return self.SS, self.JR
            

        
    def undiagonalize(self,with_remainder=True):
        """Transformed the Hamiltonian to the basis before diagonalization
        
        Parameters
        ----------
        with_remainder : bool
            Specifies if we add the coupling smaller than the cutt-off
            used in diagonalization back to the Hamiltonian.
        
        """
        self.data = numpy.dot(self.SS,numpy.dot(self.data,self.SS.T))
        if self._has_remainder_coupling:
            self.JR = numpy.dot(self.SS,numpy.dot(self.JR,self.SS.T))
        if with_remainder and self._has_remainder_coupling:                
            self.data += self.JR
            
    def remove_cutoff_coupling(self, coupling_cutoff):
        """Removes the couplings smaller than a specified cutoff
        
        Parameters
        ----------
        
        coupling_cutoff : float, optional
            Specifies the smallest (absolute) value of coupling 
            which is taken into account. Smaller couplings are removed
            and a remainder coupling matrix is returned together with
            the diagonalization matrix (see Returns section).
            
        """
        if coupling_cutoff is None:
            coupling_cutoff = 0.0
        if coupling_cutoff < 0.0:
            raise Exception("Coupling cutoff value must be positive")
            
        JR = numpy.zeros((self.dim,self.dim),dtype=numpy.float64)
        # go through all couplings and remove small ones
        for ii in range(self.dim):
            for jj in range(ii+1,self.dim):
                if (numpy.abs(self.data[ii,jj])
                        < numpy.abs(coupling_cutoff)):
                    JR[ii,jj] = self._data[ii,jj]
                    JR[jj,ii] = self._data[jj,ii]
                    self._data[ii,jj] = 0.0
                    self._data[jj,ii] = 0.0
        self.JR = JR
        self._has_remainder_coupling = True
        
    def subtract_cutoff_coupling(self, coupling_cutoff):
        """Subtracts the cut-off coupling from all coupling elements
        
        We supress the couplings by a given amount. If the coupling is
        smaller than the cutoff it is removed, if it is larger than 
        cutoff, the cutoff values is subtracted from the absolute value,
        and the cutoff values is stored in the object for subsequent
        restoration.
        
        Parameters
        ----------
        
        coupling_cutoff : float, optional
            Specifies the smallest (absolute) value of coupling 
            which is taken into account. Smaller couplings are removed
            and a remainder coupling matrix is returned together with
            the diagonalization matrix (see Returns section).
            
        """        
        if coupling_cutoff is None:
            coupling_cutoff = 0.0
            
        if type(coupling_cutoff) in (list, tuple, numpy.ndarray):
            
            raise Exception("Variable coupling cutoff not implemented yet")
        
        else:
            
            if coupling_cutoff < 0.0:
                raise Exception("Coupling cutoff value must be positive")
                
            coupcut = self.convert_2_internal_u(coupling_cutoff)
                
            JR = numpy.zeros((self.dim,self.dim),dtype=numpy.float64)
            # go through all couplings and remove small ones
            for ii in range(self.dim):
                for jj in range(ii+1,self.dim):
                    #
                    # if the coupling <= coupling_cutoff -> remove it
                    #
                    if (numpy.abs(self._data[ii,jj])
                            <= numpy.abs(coupcut)):
                        JR[ii,jj] = self._data[ii,jj]
                        JR[jj,ii] = self._data[jj,ii]
                        self._data[ii,jj] = 0.0
                        self._data[jj,ii] = 0.0
                    #
                    # if the coupling > coupling_cutoff -> suppress it
                    #
                    else:
                        absv = numpy.abs(self._data[ii,jj]) 
                        sign = self._data[ii,jj]/absv
                        val = absv - coupcut
                        if val < 0.0:
                            val = 0.0
                        JR[ii,jj] = sign*coupcut
                        self._data[ii,jj] = sign*val
                        JR[jj,ii] = JR[ii,jj]
                        self._data[jj,ii] = self._data[ii,jj]

        self.JR = JR
        self._has_remainder_coupling = True
        
        
    def recover_cutoff_coupling(self):
        """
        
        """
        if self._has_remainder_coupling:
            self._data += self.JR 
            self.JR[:,:] = 0.0
            self._has_remainder_coupling = False

    def transform(self,SS,inv=None):
        """Transformation of the Hamiltonian and the remainder coupling
        
        
        This function transforms the Operator into a different basis, using
        a given transformation matrix.
        
        Parameters
        ----------
         
        SS : matrix, numpy.ndarray
            transformation matrix
            
        inv : matrix, numpy.ndarray
            inverse of the transformation matrix
            
        """        
        if (self.manager.warn_about_basis_change):
                print("\nQr >>> Operator '%s' changes basis" %self.name)
        
        if inv is None:
            S1 = numpy.linalg.inv(SS)
        else:
            S1 = inv

        self._data = numpy.dot(S1,numpy.dot(self._data,SS))
        
        if self._has_remainder_coupling:
            self.JR = numpy.dot(S1,numpy.dot(self.JR,SS))
            
            
    def __str__(self):
        out  = "\nquantarhei.Hamiltonian object"
        out += "\n============================="
        out += "\nunits of energy %s" % self.unit_repr()
        out += "\ndata = \n"
        out += str(self.data)
        return out