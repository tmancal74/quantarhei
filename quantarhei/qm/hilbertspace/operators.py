# -*- coding: utf-8 -*-

import numbers
import numpy

from ...core.matrixdata import MatrixData
from ...core.saveable import Saveable
from ...utils.types import BasisManagedComplexArray
from ...core.managers import BasisManaged
from .statevector import StateVector
from ... import COMPLEX, REAL



class Operator(MatrixData, BasisManaged, Saveable):
    """Class representing quantum mechanical operators
    
    
    The class Operator represents an operator on quantum mechanical Hilbert
    space. Together with StateVector they provide the basic functionality 
    for quantum mechanical calculations in Hilbert space.    
    
           
    """
    
    data = BasisManagedComplexArray("data")    
    
    def __init__(self, dim=None, data=None, real=False, name=""):

        if not ((dim is None) and (data is None)):
            # Set the currently used basis
            cb = self.manager.get_current_basis()
            self.set_current_basis(cb)
            # unless it is the basis outside any context
            if cb != 0:
                self.manager.register_with_basis(cb, self)
                
            self.name=name
                 
            # set data
            if (dim is None) and (data is None):
                raise Exception() #HilbertSpaceException
            
            if data is not None:
                if isinstance(data,list):
                    data = numpy.array(data)
                if Operator.assert_square_matrix(data):
                    if dim is not None:
                        # shape 1 is OK even for time dependent operators
                        # and time dependent StateVector
                        if data.shape[1] != dim:
                            raise Exception() #HilbertSpaceException
                    self.data = data
                    self.dim = self._data.shape[1]
                else:
                    raise Exception() #HilbertSpaceException
            else:
                if real:
                    self.data = numpy.zeros((dim,dim),dtype=REAL)
                else:
                    self.data = numpy.zeros((dim,dim),dtype=COMPLEX)
                self.dim = dim


    def __add__(self, other):
        """Addition of two operators. Returns self.
        
        """
        self.data += other.data            
        return self


    def apply(self, obj):
        """Apply the operator to vector or operator on the right
        
        """
        
        if isinstance(obj, Operator):
            
            return Operator(data=numpy.dot(self.data,obj.data))
        
        elif isinstance(obj, StateVector):
        
            return StateVector(data=numpy.dot(self.data,obj.data))
            
        else:
            
            raise Exception("Cannot apply operator to the object")
        

    def transform(self, SS, inv=None):
        """Transformation of the operator by a given matrix
        
        
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
 
    def is_diagonal(self):
        dat = self._data.copy()
        for i in range(self.dim):
            dat[i,i] = 0.0
        return numpy.allclose(dat, numpy.zeros(self.dim))
        
 
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
    
    def __init__(self, dim=None, data=None, name=""):
        
        if not ((dim is None) and (data is None)):
            Operator.__init__(self,dim=dim,data=data,name=name)
            if not self.check_selfadjoint():
                raise Exception("The data of this operator have"+
                "to be represented by a selfadjoint matrix") 
        
        
    def check_selfadjoint(self):
        return (numpy.isclose(numpy.transpose(numpy.conj(self.data)),
                              self.data)).all()
        
    def diagonalize(self):
        # first use is of "data", the rest of "_data"
        dd,SS = numpy.linalg.eigh(self.data)
        self._data = numpy.zeros(self._data.shape)
        for ii in range(self._data.shape[0]):
            self._data[ii,ii] = dd[ii]
        return SS
        
    def get_diagonalization_matrix(self):
        dd, SS = numpy.linalg.eigh(self._data)
        return SS        
    
    def __str__(self):
        out  = "\nquantarhei.SelfAdjointOperator object"
        out += "\n====================================="
        out += "\ndata = \n"
        out += str(self.data)
        return out        
        
   
class UnityOperator(SelfAdjointOperator):
    
    def __init__(self, dim=None, name=""):
        if dim is None:
            raise Exception("Dimension parameters 'dim' has to be specified")
        series = numpy.array([1 for i in range(dim)], dtype=COMPLEX)
        data = numpy.diag(series)
        super().__init__(dim=dim,data=data,name=name)    
    

class BasisReferenceOperator(SelfAdjointOperator):
    
    def __init__(self, dim=None, name=""):
        if dim is None:
            raise Exception("Dimension parameters 'dim' has to be specified")
        series = numpy.array([i for i in range(dim)], dtype=numpy.float64)
        data = numpy.diag(series)
        super().__init__(dim=dim,data=data,name=name)
    
     
class ProjectionOperator(Operator):
    """
    Projection operator projecting from state m to state n.
    Braket definition:
    
    |n\rangle \langle m|

    """    
    def __init__(self, to_state=-1, from_state=-1, dim=0):
        # here the operator is create and it will know about basis
                    
        if dim > 0:
            super().__init__(dim=dim, real=True)
            if ((to_state >= 0) and (to_state < dim)) \
                and ((from_state >= 0) and (from_state < dim)):
                self.data[to_state, from_state] = 1.0
            #else:
            #    raise Exception("Indices out of range")
        else:
            raise Exception("Wrong operator dimension")
            
            
    def __mult__(self, other):
        """Multiplication of operator by scalar
        
        """
        if isinstance(other, numbers.Number):
            self._data = self._data*other
        else:
            raise Exception("Only multiplication by scalar is allowed")
        return self._data

    def __rmult__(self, other):
        """Multiplication from right
        
        """
        return self.__mult__(other)

        
            
class DensityMatrix(SelfAdjointOperator, Saveable):
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


    def normalize(self):
        """Normalize the trace of the density matrix
        
        """
        tr = numpy.trace(self.data)
        self.data = self.data/tr


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