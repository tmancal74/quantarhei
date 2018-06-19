# -*- coding: utf-8 -*-
"""Class representing secular Superoperators



    Class Details
    -------------

"""

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
    
    >>> A.secularize()
    
    >>> A.is_secular
    True

    Nothing will be returned below (returns None)
    
    >>> A.get_secular_basis_operator()

    If we secularize in certain basis, we remember the operator which 
    defines the basis

    >>> O1 = qr.qm.SelfAdjointOperator(data=[[0.0, 1.0],[1.0, 0.0]])
    >>> B = Test()
    >>> with qr.eigenbasis_of(O1):
    ...     B.secularize()
    
    >>> B.is_secular
    True
    
    >>> O2 = B.get_secular_basis_operator()
    >>> O1 == O2
    True
    
    """    
    
    is_secular = False
    secular_basis_op = None

    ###########################################################################
    #
    #  Public methods
    #
    ###########################################################################
        
    def secularize(self):
        """Secularizes the SuperOperator in current basis
        
        
        
        
        """
        
        if not self.is_secular:
            
            self.secular_basis_op = self._get_current_basis_op()
        
            if self.as_operators:
                self._set_population_rates_from_operators()
                self._set_dephasing_rates_from_operators()        
            else:
                self._set_population_rates_from_tensor()
                self._set_dephasing_rates_from_tensor()
            
            self.is_secular = True


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
    
    
    ###########################################################################
    #
    # Private methods
    #
    ###########################################################################
    
    def _get_current_basis_op(self):
        """Returns the operator which defines the currently used basis
        
        """
        return Manager().current_basis_operator
    
    

    
    
