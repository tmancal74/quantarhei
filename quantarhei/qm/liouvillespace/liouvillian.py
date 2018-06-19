# -*- coding: utf-8 -*-
"""Liouvillian super operator

    Louvillean super operator represent the action of a commutator with
    Hamiltonian
    
    Class Details
    -------------
    
"""
from .superoperator import SuperOperator
from ..hilbertspace.operators import UnityOperator

class Liouvillian(SuperOperator):
    """Class representing and creating Liouvillian super operator
    
    
    Parameters
    ----------
    
    ham : Hamiltonian
        Hamiltonian from which we create the Liouvillian
    
    
    Examples
    --------
    
    >>> import numpy
    >>> import quantarhei as qr
    >>> HH = qr.qm.Hamiltonian(data=[[0.0, 1.0], [1.0, 0.2]])
    >>> rho = qr.qm.DensityMatrix(data=[[0.9, 0.05], [0.05, 0.1]])
    >>> Li = Liouvillian(HH)
    >>> com1 = numpy.dot(HH.data, rho.data) - numpy.dot(rho.data, HH.data)
    >>> com2 = Li.apply(rho)
    >>> numpy.allclose(com1, com2.data)
    True
    
    """
    
    def __init__(self, ham):
        super().__init__(dim=ham.dim)
        dim = ham.dim
        delta = UnityOperator(dim=dim)
        for k_i in range(dim):
            for k_j in range(dim):
                for k_k in range(dim):
                    for k_l in range(dim):
                        self.data[k_i, k_j, k_k, k_l] = \
                        ham.data[k_i, k_k]*delta.data[k_j, k_l] - \
                        delta.data[k_i, k_k]*ham.data[k_l, k_j]
        
        
