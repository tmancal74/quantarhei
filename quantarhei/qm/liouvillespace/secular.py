# -*- coding: utf-8 -*-
"""Class representing secular Superoperators



    Class Details
    -------------

"""
import numpy

from ...core.managers import Manager

class Secular:
    """Class representing secular Superoperators


    
    
    Attributes
    ----------
    
    is_secular : bool
        If True, the superoperator has a secular form and it should be
        handled like that.
        
        
    Methods
    -------
    
    secularize() :
        Converts the superoperator into secular form
   
    
    Examples
    --------
    
    This class is used as a base class for classis that want to provide 
    a secular form. The class has to define the attribute as_operators (bool)
    and it has to be BasisManaged
    
    >>> import quantarhei as qr
    >>> class Test(Secular, qr.core.BasisManaged):
    ...     def __init__(self):
    ...         self.as_operators = False
    ...         self.data = None
    ...     def _set_population_rates_from_tensor(self):
    ...         pass
    ...     def _set_dephasing_rates_from_tensor(self):
    ...         pass
    
    >>> A = Test()
    >>> A.is_secular
    False
    
    >>> A.secularize(use_data=False)
    
    >>> A.is_secular
    True

    Nothing will be returned below (returns None)
    
    >>> A.get_secular_basis_operator()

    If we secularize in certain basis, we remember the operator which 
    defines the basis

    >>> O1 = qr.qm.SelfAdjointOperator(data=[[0.0, 1.0],[1.0, 0.0]])
    >>> B = Test()
    >>> with qr.eigenbasis_of(O1):
    ...     B.secularize(use_data=False)
    
    >>> B.is_secular
    True
    
    >>> O2 = B.get_secular_basis_operator()
    >>> O1 == O2
    True
    
    """    
    
    # True if the object is secular
    is_secular = False
    
    # True if the data are time-dependent
    secular_time_dependent = False
    
    # True if data data property is secular
    secular_data = False
    
    # True if the data can be used to remove secularization
    secular_reversible = False
    
    # operator whose eigenstates for the basis in which this object is secular
    secular_basis_op = None
    

    #
    # Data are stored in these two properties if the object is secular
    #
    
    # coherence decay rates
    secular_GG = None
    
    # population tranfer rates
    secular_KK = None

    ###########################################################################
    #
    #  Public methods
    #
    ###########################################################################
        
    def secularize(self, reversible=False, use_data=True):
        """Secularizes the SuperOperator in current basis
        
        
        Parameters
        ----------
        
        reversible : bool
            If the secularization is reversible, the original data property
            is left untouched, so that secularization can be reversed.
        
        
        """
        if reversible and use_data:
            raise Exception("data cannot be secularized reversibly")
            
        self.secular_reversible = reversible
        
        if use_data:
            if self.data is None:
                raise Exception("Cannot use data when data is None")
            self.secular_data = True
        
        if not self.is_secular:
            
            self.secular_basis_op = self._get_current_basis_op()
            
            if self.secular_data:
                
                # do it to data
                self._secularize_data()
            
            else:
                
                # create rate matrices
                if self.as_operators:
                    self._set_population_rates_from_operators()
                    self._set_dephasing_rates_from_operators()        
                else:
                    self._set_population_rates_from_tensor()
                    self._set_dephasing_rates_from_tensor()
            
            self.is_secular = True


    def recover_nonsecular(self):
        
        if self.is_secular:
            if self.secular_reversible:
            
                self.is_secular = False
                self.secular_basis_op = None
                self.secular_reversible = False
                self.secular_data = False
                self.secular_time_dependent = False
                self.secular_GG = None
                self.secular_KK = None
            
            else:
                
                raise Exception("Cannot recover non-secular data:"
                                +" secularization was performed irreversibly")
             
            
            

    def get_secular_basis_operator(self):
        """Returns the operator to ensure the correct basis context
        
        """
        return self.secular_basis_op


    ###########################################################################
    #
    # These methods need to be implemented in classes which inherit from
    # Secular. 
    #
    ###########################################################################

    def _set_population_rates_from_operators(self):
        raise Exception("Method _set_population_rates_from_operators()"+
                        " needs to be implemented in a class inheriting"+
                        " from Secular")


    def _set_population_rates_from_tensor(self):
        raise Exception("Method _set_population_rates_from_tensor()"+
                        " needs to be implemented in a class inheriting"+
                        " from Secular")


    def _set_dephasing_rates_from_operators(self):
        raise Exception("Method _set_dephasing_rates_from_operators()"+
                        " needs to be implemented in a class inheriting"+
                        " from Secular")


    def _set_dephasing_rates_from_tensor(self):
        raise Exception("Method _set_dephasing_rates_from_tensor()"+
                        " needs to be implemented in a class inheriting"+
                        " from Secular")


    def _secularize_data(self):

        # FIXME: temporary fix
        if self.as_operators:
            self.convert_2_tensor()   
            
        if self.data.ndim == 4:
            N = self.data.shape[0]
            for ii in range(N):
                for jj in range(N):
                    for kk in range(N):
                        for ll in range(N):
                            if not (((ii == jj) and (kk == ll)) 
                                or ((ii == kk) and (jj == ll))) :
                                    self.data[ii,jj,kk,ll] = 0
        else:  
            N = self.data.shape[1]
            for ii in range(N):
                for jj in range(N):
                    for kk in range(N):
                        for ll in range(N):
                            if not (((ii == jj) and (kk == ll)) 
                                or ((ii == kk) and (jj == ll))) :
                                    self.data[:,ii,jj,kk,ll] = 0
        
    
    def apply(self, rho):
        """Application of the secular tensor on the statistical operator
        
        
        """
#        rhoret = -self.secular_GG*rho.data
#        if self._has_rates:
#            rhoret += numpy.diagflat(numpy.einsum("ij,jj", 
#                                     self.secular_KK.data, rho.data))
            
        # FIXME: time this to find out how bad is to call additional function
        return self._apply(rho.data)
        

    def _apply(self, rho):
        """Application of the tensor directly to an array
        
        """
        rhoret = -self.secular_GG*rho
        if self._has_rates:
            rhoret += numpy.diagflat(numpy.einsum("ij,jj", 
                                     self.secular_KK.data, rho))
        
        
    ###########################################################################
    #
    # Private methods
    #
    ###########################################################################
    
    def _get_current_basis_op(self):
        """Returns the operator which defines the currently used basis
        
        """
        return Manager().current_basis_operator
    
    

    
    
