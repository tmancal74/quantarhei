# -*- coding: utf-8 -*-

from .operators import SelfAdjointOperator
from ...core.managers import BasisManaged
from ...utils.types import BasisManagedReal

import numpy


class Hamiltonian(SelfAdjointOperator,BasisManaged):
    """Hamiltonian operator
    
    
    
    """
    
    _has_remainder_coupling = False
    
    data = BasisManagedReal("data")
    
    def diagonalize(self,coupling_cutoff=None):
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
            return SS
        else:
            JR = numpy.zeros((self.dim,self.dim),dtype=numpy.float64)
            # go through all couplings and remove small ones
            for ii in range(self.dim):
                for jj in range(ii+1,self.dim):
                    if (numpy.abs(self.data[ii,jj])
                        < numpy.abs(coupling_cutoff)):
                            JR[ii,jj] = self.data[ii,jj]
                            JR[jj,ii] = JR[ii,jj]
                            self.data[ii,jj] = 0.0
                            self.data[jj,ii] = 0.0
            # diagonalize the strong coupling part
            dd,SS = numpy.linalg.eigh(self.data)
            self.data = numpy.zeros(self.data.shape,dtype=numpy.float64)
            for ii in range(0,self.data.shape[0]):
                self.data[ii,ii] = dd[ii]
            # transform the remainder of couling correspondingly
            JR = numpy.dot(SS.T,numpy.dot(JR,SS))
            self.SS = SS
            self.JR = JR
            self._has_remainder_coupling = True
            return SS,JR
            

        
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
        if with_remainder:                
            self.data += self.JR
            
            
    def __str__(self):
        out  = "\nquantarhei.Hamiltonian object"
        out += "\n============================="
        out += "\ndata = \n"
        out += str(self.data)
        return out