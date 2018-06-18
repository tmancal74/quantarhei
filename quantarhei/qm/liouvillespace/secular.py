# -*- coding: utf-8 -*-
"""Class representing secular Superoperators



    Class Details
    -------------

"""


class Secular:
    """Class representing secular Superoperators

    

   
    
    """    
    
    is_secular = False

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
        
            if self._as_operators:
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
        
        return None
    
    

    
    
