# -*- coding: utf-8 -*-

from ...core.matrixdata import MatrixData
import numpy
import scipy


class Operator(MatrixData):
    """
    The class Operator represents an operator on quantum mechanical Hilbert space.
    Together with StateVector they provide the basic functionality for quantum
    mechanical calculations in Hilbert space.    
    
    The Operator provides two constants, values 0 and 1, to identify
    the storage mode of the data. The following calls reveal their values.
    
    >>> Operator.STORAGE_MODE_VECTOR
    0
    >>> Operator.STORAGE_MODE_MATRIX
    1
    
    Storage mode is by default set to matrix
    
    >>> o = Operator(2)
    >>> o.getStorageMode()
    1
    >>> o.data
    array([[ 0.,  0.],
           [ 0.,  0.]])

    >>> o.data[0,0] = 0.0
    >>> o.data[1,0] = 1.0
    >>> o.data[0,1] = 0.1
    >>> o.data[1,1] = 1.1
    >>> o.setStorageMode(Operator.STORAGE_MODE_VECTOR)
    >>> o.data
    array([ 0. ,  0.1,  1. ,  1.1])
        
    >>> o.setStorageMode(Operator.STORAGE_MODE_MATRIX)
    >>> o.data
    array([[ 0. ,  0.1],
           [ 1. ,  1.1]])
           
    """
    STORAGE_MODE_VECTOR = 0
    STORAGE_MODE_MATRIX = 1
    
    def __init__(self,dim=None,data=None,
                 storage_mode=STORAGE_MODE_MATRIX,real=False):
        """
        
        Operator can be initiallized either by data or by setting the dimension.
        In the latter case, its elements will be set zero. Alternatively,
        storage mode can be specified.
        
        >>> o = Operator(dim=3)
        >>> o.data.shape
        (3, 3)
        
        or 
        
        >>> A = numpy.array([[1.0, 0.1],[0.1, 1.0]])
        >>> o = Operator(data=A)
        >>> o.data.shape
        (2, 2)

        Operator must be a square matrix of numpy.ndarray type. Setting
        an non-square matrix or a list as imput data will lead to an Exception
    
        >>> A = [[1.0, 0.0]]
        >>> o = Operator(data=A)
        Traceback (most recent call last):
            ...
            raise HilbertSpaceException
        cu.oqs.hilbertspace.HilbertSpaceException
        
        The problem is not elevated even if you make it square
        
        >>> A = [[1.0, 0.0],[0.0, 0.0]]
        >>> o = Operator(data=A)
        Traceback (most recent call last):
            ...
            raise HilbertSpaceException
        cu.oqs.hilbertspace.HilbertSpaceException

        Here is another failure:
        
        >>> A = numpy.array([[1.0, 0.0],[0.0, 0.0]])
        >>> o = Operator(3,A)
        Traceback (most recent call last):
            ...
            raise HilbertSpaceException
        cu.oqs.hilbertspace.HilbertSpaceException

        
        """

        if (dim is None) and (data is None):
            raise Exception() #HilbertSpaceException
        
        if data is not None:
            if Operator.assert_square_matrix(data):
                if dim is not None:
                    if data.shape[0] != dim:
                        raise Exception() #HilbertSpaceException
                self.data = data
                self.dim = self.data.shape[0]
            else:
                raise Exception #HilbertSpaceException
        else:
            if real:
                self.data = numpy.zeros((dim,dim),dtype=numpy.float64)
            else:
                self.data = numpy.zeros((dim,dim),dtype=numpy.complex128)
            self.dim = dim
            
        self.storage_mode = Operator.STORAGE_MODE_MATRIX 
        
        self.setStorageMode(storage_mode)
        



    def transform(self,SS):
        """
        This function transforms the Operator into a different basis, using
        a given transformation matrix.
        """
        self.setStorageMode(Operator.STORAGE_MODE_MATRIX)
        S1 = scipy.linalg.inv(SS)
        self.data = numpy.dot(S1,numpy.dot(self.data,SS))
        
    def getStorageMode(self):
        return self.storage_mode
        

    def setStorageMode(self,mode):
        
        if self.storage_mode != mode:
            
            if mode == self.STORAGE_MODE_VECTOR:

                sz = self.data.shape[0]*self.data.shape[1]
                self.data.shape = (sz)
                self.storage_mode = mode
                
            if mode == self.STORAGE_MODE_MATRIX:
                
                sz = numpy.sqrt(self.data.shape[0])
                self.data.shape = (sz,sz)
                self.storage_mode = mode
                
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
        return str(self.data)
         
           
class SelfAdjointOperator(Operator):
    """
    The class SelfAdjointOperator extends the Operaator class by automatically
    checking the data for the self adjoint property. Physically measurable
    quantities are represented by self adjoint operators in quantum mechanics.
    
    >>> dm = SelfAdjointOperator(dim=3)
    >>> dm.data.shape
    (3, 3)
    
    >>> A = numpy.array([[1.0, 0.1],[0.1, 1.0]])
    >>> o = SelfAdjointOperator(data=A)
    >>> o.data.shape
    (2, 2)
        
    """
    
    def __init__(self,dim=None,data=None,storage_mode=Operator.STORAGE_MODE_MATRIX):
        Operator.__init__(self,dim=dim,data=data,storage_mode=storage_mode)
        if not self.check_selfadjoint():
            raise Exception("The data of this operator have to be represented by a selfadjoint matrix") 
        
        
    def check_selfadjoint(self):
        return (numpy.isclose(numpy.transpose(numpy.conj(self.data)),
                              self.data)).all()
        
    def diagonalize(self):
        dd,SS = numpy.linalg.eigh(self.data)
        self.data = numpy.zeros(self.data.shape)
        for ii in range(0,self.data.shape[0]):
            self.data[ii,ii] = dd[ii]
        return SS
        
        
class ProjectionOperator(Operator):
    """
    Projection operator projecting from state m to state n.
    Braket definition:
    
    |n\rangle \langle m|

    """    
    def __init__(self,n,m,dim=0,storage_mode=Operator.STORAGE_MODE_MATRIX):
        Operator.__init__(self,dim=dim,data="",storage_mode=storage_mode)
        if (n < dim and m < dim):
            self.data[n,m] = 1
        else:
            raise Exception("Projection Operator indices exceed its dimension")
            
class DensityMatrix(SelfAdjointOperator):
    """
    Density matrix is a good example of self adjoint operator. It extends
    the class SelfAdjointOperator without adding any functionality.    
    
    >>> dm = DensityMatrix(dim=3)
    >>> dm.data.shape
    (3, 3)
    
    >>> A = numpy.array([[1.0, 0.1],[0.1, 1.0]])
    >>> o = DensityMatrix(data=A)
    >>> o.data.shape
    (2, 2)
        
    """
    
    def normalize2(self,norm=1.0):
        tr = numpy.trace(self.data)
        self.data = self.data*(norm/tr)