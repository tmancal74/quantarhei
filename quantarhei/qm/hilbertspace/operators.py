# -*- coding: utf-8 -*-

from ...core.matrixdata import MatrixData
from ...utils.types import BasisManagedComplex
from ...core.managers import BasisManaged


import numpy
#import scipy


class Operator(MatrixData,BasisManaged):
    """Class representing quantum mechanical operators
    
    
    The class Operator represents an operator on quantum mechanical Hilbert
    space. Together with StateVector they provide the basic functionality 
    for quantum mechanical calculations in Hilbert space.    
    
           
    """
    
    data = BasisManagedComplex("data")    
    
    def __init__(self,dim=None,data=None,real=False):

        # Set the currently used basis
        cb = self.manager.get_current_basis()
        self.set_current_basis(cb)
        # unless it is the basis outside any context
        if cb != 0:
            self.manager.register_with_basis(cb,self)
             
        # set data
        if (dim is None) and (data is None):
            raise Exception() #HilbertSpaceException
        
        if data is not None:
            if isinstance(data,list):
                data = numpy.array(data)
            if Operator.assert_square_matrix(data):
                if dim is not None:
                    if data.shape[0] != dim:
                        raise Exception() #HilbertSpaceException
                self._data = data
                self.dim = self._data.shape[0]
            else:
                raise Exception #HilbertSpaceException
        else:
            if real:
                self._data = numpy.zeros((dim,dim),dtype=numpy.float64)
            else:
                self._data = numpy.zeros((dim,dim),dtype=numpy.complex128)
            self.dim = dim



    def transform(self,SS,inv=None):
        """
        This function transforms the Operator into a different basis, using
        a given transformation matrix.
        """        
        if inv is None:
            S1 = numpy.linalg.inv(SS)
        else:
            S1 = inv

        #S1 = scipy.linalg.inv(SS)
        self._data = numpy.dot(S1,numpy.dot(self._data,SS))
        
                
    def assert_square_matrix(A):            
        if isinstance(A,numpy.ndarray):
            if A.ndim == 2:
                if A.shape[0] == A.shape[1]:
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False
 
    def __str__(self):
        out  = "\nquantarhei.Operator object"
        out += "\n=========================="
        out += "\ndata = \n"
        out += str(self.data)
        return out
         
           
class SelfAdjointOperator(Operator):
    """Class representing a self adjoint operator
    
    The class SelfAdjointOperator extends the Operator class by automatically
    checking the data for the self adjoint property. Physically measurable
    quantities are represented by self adjoint operators in quantum mechanics.
    
        
    """
    
    def __init__(self,dim=None,data=None):
        Operator.__init__(self,dim=dim,data=data)
        if not self.check_selfadjoint():
            raise Exception("The data of this operator have"+
            "to be represented by a selfadjoint matrix") 
        
        
    def check_selfadjoint(self):
        return (numpy.isclose(numpy.transpose(numpy.conj(self.data)),
                              self.data)).all()
        
    def diagonalize(self):
        dd,SS = numpy.linalg.eigh(self.data)
        self.data = numpy.zeros(self.data.shape)
        for ii in range(0,self.data.shape[0]):
            self.data[ii,ii] = dd[ii]
        return SS
        
    def __str__(self):
        out  = "\nquantarhei.SelfAdjointOperator object"
        out += "\n====================================="
        out += "\ndata = \n"
        out += str(self.data)
        return out        
        
        
        
class ProjectionOperator(Operator):
    """
    Projection operator projecting from state m to state n.
    Braket definition:
    
    |n\rangle \langle m|

    """    
    def __init__(self,n,m,dim=0):
        Operator.__init__(self,dim=dim,data="")
        if (n < dim and m < dim):
            self.data[n,m] = 1
        else:
            raise Exception("Projection Operator indices exceed its dimension")
            
            
            
            
class DensityMatrix(SelfAdjointOperator):
    """Class representing a density matrix
    
    
    Density matrix is a good example of self adjoint operator. It extends
    the class SelfAdjointOperator without adding much functionality.    

        
    """
    
    
    def normalize2(self,norm=1.0):
        # using self.data to allow transformation
        tr = numpy.trace(self.data)
        # using data because any transformation needed was already done
        self._data = self._data*(norm/tr)
        
    def get_populations(self):
        """Returns a vector of diagonal elements of the density matrix
        
        """
        # using self.data to allow transformation of the basis
        pvec = numpy.zeros(self.data.shape[0],dtype=numpy.float)
        # using self._data because transformations are irrelevant 
        for n in range(self._data.shape[0]):
            # using self._data because any transformation needed was already
            # made above
            pvec[n] = numpy.real(self._data[n,n])
            
        return pvec
        
    def excite_delta(self,dmoment,epolarization=[1.0, 0.0, 0.0]):
        """Returns a density matrix obtained by delta-pulse excitation
        
        """
        # using .data to allow transformation 
        dd = dmoment.data
        etimesd = numpy.zeros((dd.shape[0],dd.shape[1]),dtype=numpy.float)
        for i in range(3):
            etimesd += dd[:,:,i]*epolarization[i]
            
        # using self.data to allow transformation to current basis
        dat = numpy.dot(etimesd,numpy.dot(self.data,etimesd))
        
        return ReducedDensityMatrix(data=dat)
        
            
    def __str__(self):
        out  = "\nquantarhei.DensityMatrix object"
        out += "\n==============================="
        out += "\ndata = \n"
        out += str(self.data)
        return out             
        
        
        
class ReducedDensityMatrix(DensityMatrix):
    """Class representing a reduced density matrix
    
    This class is basically just a nickname for the Density Matrix
    
    
    """
    
    def __str__(self):
        out  = "\nquantarhei.ReducedDensityMatrix object"
        out += "\n======================================"
        out += "\ndata = \n"
        out += str(self.data)
        return out  