 # -*- coding: utf-8 -*-
"""Class adding exciton analysis on molecular aggregates

    Most of the present methods are available after the aggregate is
    diagonalized by calling the ``diagonalize`` method.

    **This class should not be used directly**. Use `Aggregate` class, which
    inherits all the methods from here, instead.
    
    Examples
    --------
    
    >>> import quantarhei as qr
    >>> agg = qr.TestAggregate("homodimer-2")
    >>> agg.set_coupling_by_dipole_dipole()
    >>> agg.build()
    >>> # create information about eigenstates of the aggregate (to call `diagonalize()` is crucial)
    >>> agg.diagonalize()
    >>> #
    >>> # Create a report on expantions of state with index `1`
    >>> agg.report_on_expansion(state=1)
    +-------+-------------+--------------+------------------+
    | index | squares     | coefficients | state signatures |
    +-------+-------------+--------------+------------------+
    | 1     | 0.50000000  | -0.70710678  | ((1, 0), ())     |
    | 2     | 0.50000000  | 0.70710678   | ((0, 1), ())     |
    | 0     | 0.00000000  | 0.00000000   | ((0, 0), ())     |
    | 0     | -1.00000000 | 0.00000000   | ((0, 0), ())     |
    | 0     | -1.00000000 | 0.00000000   | ((0, 0), ())     |
    +-------+-------------+--------------+------------------+
    

    Class Details
    -------------

"""
import numpy

from .aggregate_spectroscopy import AggregateSpectroscopy
import quantarhei as qr

class AggregateExcitonAnalysis(AggregateSpectroscopy):
    """Class adding exciton analysis on molecular aggregates
    
    
    
    """
    
    
    def get_expansion_squares(self, state=0):
        """Returns the squares of expansion coefficients of an excitonic state.

        The Aggregate must be built and diagonalized. This output is used by
        the `report_on_expansion()` method.
        
        >>> import quantarhei as qr
        >>> agg = qr.TestAggregate("homodimer-2")
        >>> agg.set_coupling_by_dipole_dipole()
        ...
        ... # build and diagonalize
        >>> agg.build()  
        >>> agg.diagonalize()
        ...
        ... # expansion coefficients 
        >>> (indx, coefs, sqrs) = agg.get_expansion_squares(1) 
        >>> print(indx)
        [0, 1, 2]
        >>> print(coefs)
        [ 0.         -0.70710678  0.70710678]
        >>> print(sqrs)
        [ 0.   0.5  0.5]
        
        
        """
        coefs = self.SS[:,state]
        sqrs  = self.SS[:,state]**2
        indx = [i for i in range(self.HH.shape[0])]
        return indx, coefs, sqrs
        
        
    def report_on_expansion(self, file=None, state=0, N=5):
        """Prints a short report on the composition of an exciton state

        Parameters
        ----------
        
        state : int
            Excitonic state of the aggregate
            
        N : int
            Number of states in expansion to report
            
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> agg = qr.TestAggregate("dimer-2")
        >>> agg.set_coupling_by_dipole_dipole()
        ...
        ... # build and diagonalize
        >>> agg.build()  
        >>> agg.diagonalize()
        ...
        >>> agg.report_on_expansion(state=1)
        +-------+-------------+--------------+------------------+
        | index | squares     | coefficients | state signatures |
        +-------+-------------+--------------+------------------+
        | 1     | 0.90998848  | -0.95393316  | ((1, 0), ())     |
        | 2     | 0.09001152  | 0.30001920   | ((0, 1), ())     |
        | 0     | 0.00000000  | 0.00000000   | ((0, 0), ())     |
        | 0     | -1.00000000 | 0.00000000   | ((0, 0), ())     |
        | 0     | -1.00000000 | 0.00000000   | ((0, 0), ())     |
        +-------+-------------+--------------+------------------+
    
        """
        try:
            from terminaltables import AsciiTable
        except:
            raise Exception("Get terminaltables package "
                            +"for this functionality")
    
        indx, coefs, sqrs = self.get_expansion_squares(state)
    
        table_data = []
        table_data.append(["index","squares", "coefficients",
                           "state signatures"])
        for i in range(N):
            imax, sqr = _strip_max_coef(indx, sqrs)
            coef = coefs[imax]

            sqr_s = "{0:.8f}".format(sqr)
            coef_s = "{0:.8f}".format(coef)
            
            sta = self.get_state_signature_by_index(imax)
            table_data.append([imax, sqr_s, coef_s, sta])
            
            #print(imax, "\t", coef,"\t", sta)
            
        table = AsciiTable(table_data)
        print(table.table, file=file)


    def get_intersite_mixing(self, state1=0, state2=0):
        """Returns inter site mixing ration
        
        Inter site mixing ratio gives the probability that
        states ``state1`` and ``state2`` are localized
        on different pigments. This measure can be used to distinguish 
        different types of coherence, e.g. vibrational from purely electronic.
        
        **Literature:** 
        
        |Cit1| available at |citlink|_
        
        .. |Cit1| replace:: Pavel Maly, Oscar J. G. Somsen, Vladimir I. \
        Novoderezhkin, Tomas Mancal, and Rienk van Grondelle, *The Role \
        of Resonant Vibrations in Electronic Energy Transfer*, ChemPhysChem \
        2016, 17,1356–1368
        
        .. |citlink| replace:: DOI:10.1002/cphc.201500965 
        .. _citlink: http://dx.doi.org/10.1002/cphc.201500965
        
         
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> agg = qr.TestAggregate("dimer-2")
        >>> agg.set_coupling_by_dipole_dipole()
        ...
        ... # build and diagonalize
        >>> agg.build()  
        >>> agg.diagonalize()

        >>> print("{0:.4f}".format(agg.get_intersite_mixing(1,2)))   
        0.8362
        
        """
        ci = self.SS[:,state1]
        cj = self.SS[:,state2]
        
        xi = 1.0
        
        # loop over all electronic states
        for n_el in range(self.Nel):
            # loop over all vibrational states in a given electronic state
            for n_nu in self.vibindices[n_el]:
                for n_mu in self.vibindices[n_el]:
                    x = (ci[n_nu]**2)*(cj[n_mu]**2)
                    xi -= x
                    
        return xi    


    def get_transition_dipole(self, state1=0, state2=0):
        """Returns transition dipole moment between two state.
 
        If the second state is not specified, we get transition to the first
        state from the ground state.
        
        Parameters
        ----------
        
        state1 : int
            Starting state of the transition 
            
        state2 : int
            Final states of the transition
            
        
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> agg = qr.TestAggregate("dimer-2")
        >>> agg.set_coupling_by_dipole_dipole()
        ...
        ... # build and diagonalize
        >>> agg.build()  
        >>> agg.diagonalize() 
        ...
        >>> a = agg.get_transition_dipole(0,1)
        >>> print("{0:.6f}".format(a))
        2.480265
        
        
        """
        if self._diagonalized:
            return self.D2[state1, state2]
        
        else:
            raise Exception("Aggregate has to be diagonalized")
            
            
    def get_state_energy(self, state=0):
        """Return the energy of a state with a given index
        
        
        Parameters
        ----------
        
        state : int
            Index of the state whose energy we ask for
            
        
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> agg = qr.TestAggregate("dimer-2")
        >>> agg.set_coupling_by_dipole_dipole()
        ...
        ... # build and diagonalize
        >>> agg.build()  
        >>> agg.diagonalize()
        ...
        >>> print("{0:.4f}".format(agg.get_state_energy(2)))
        2.3231
        
        >>> with qr.energy_units("1/cm"):
        ...     print("{0:.4f}".format(agg.get_state_energy(2)))
        12332.9320
        
        
        """
        if self._diagonalized:
            return self.convert_energy_2_current_u(self.HD[state]) 
        else:
            raise Exception("Aggregate has to be diagonalized")


    def exciton_report(self, file=None, start=1, stop=None, Nrep=5, criterium=None):
        """Prints a report on excitonic properties of the aggregate
        
        Parameters
        ----------
        
        start : int (default = 1)
            First state to report on. Because 0 is often the ground state,
            it is skipped. For systems with molecular vibrations, first mixed
            states start with an index different from 1.
            
        stop : int (default = None)
            Where to stop the report. if None, all states are reported
            
        Nrep : int (default = 5)
            How many sites of the exciton decomposition we should list   
            
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> agg = qr.TestAggregate("dimer-2")
        >>> agg.set_coupling_by_dipole_dipole()
        ...
        ... # build and diagonalize
        >>> agg.build()  
        >>> agg.diagonalize()
        ...
        >>> agg.exciton_report()
        Report on excitonic properties
        ------------------------------
        <BLANKLINE>
        Exciton 1
        =========
        <BLANKLINE>
        Transition energy        : 11967.06803045 1/cm
        Transition dipole moment : 2.48026499 D
        +-------+-------------+--------------+------------------+
        | index | squares     | coefficients | state signatures |
        +-------+-------------+--------------+------------------+
        | 1     | 0.90998848  | -0.95393316  | ((1, 0), ())     |
        | 2     | 0.09001152  | 0.30001920   | ((0, 1), ())     |
        | 0     | 0.00000000  | 0.00000000   | ((0, 0), ())     |
        | 0     | -1.00000000 | 0.00000000   | ((0, 0), ())     |
        | 0     | -1.00000000 | 0.00000000   | ((0, 0), ())     |
        +-------+-------------+--------------+------------------+
        <BLANKLINE>
        Exciton 2
        =========
        <BLANKLINE>
        Transition energy        : 12332.93196955 1/cm
        Transition dipole moment : 5.16973501 D
        +-------+-------------+--------------+------------------+
        | index | squares     | coefficients | state signatures |
        +-------+-------------+--------------+------------------+
        | 2     | 0.90998848  | 0.95393316   | ((0, 1), ())     |
        | 1     | 0.09001152  | 0.30001920   | ((1, 0), ())     |
        | 0     | 0.00000000  | 0.00000000   | ((0, 0), ())     |
        | 0     | -1.00000000 | 0.00000000   | ((0, 0), ())     |
        | 0     | -1.00000000 | 0.00000000   | ((0, 0), ())     |
        +-------+-------------+--------------+------------------+
        <BLANKLINE>        
        """
        
        start_at = start
        stop_at = stop
        
        
        print("Report on excitonic properties", file=file)
        print("------------------------------\n", file=file)
        N01 = self.Nb[0]+self.Nb[1]
        if self.mult > 1:
            N01 += self.Nb[2]
            
        for Nst in range(N01):
            
            if stop_at is None:
                stop_at = N01 + 1
                
            if (Nst >= start_at) and (Nst <= stop_at):
                
                with qr.energy_units("1/cm"):
                    tre = self.get_state_energy(Nst) - self.get_state_energy(0)
                dip = self.get_transition_dipole(0, Nst)
                
                if criterium is not None:
                    cond = criterium([tre, dip])
                    
                else:
                    cond = True
                
                if cond:
                    txt = "Exciton "+str(Nst)
                    Nlength = len(txt)
                    line = "="*Nlength
                    print(txt, file=file)
                    print(line, file=file)
                    print("", file=file)
                    with qr.energy_units("1/cm"):
                        print("Transition energy        "+
                              ": {:.8f} 1/cm".format(tre), file=file)
                        print("Transition dipole moment "+
                              ": {:.8f} D".format(dip), file=file)
                    self.report_on_expansion(file=file, state=Nst, N=Nrep)
                    print("", file=file)


    #
    #  The following routine does not work with excitons, but rather with
    #  site basis representation 
    #

    # FIXME: move it somewhere else
    def get_state_signature_by_index(self, N):
        """Return aggregate vibronic state signature  by its index
        
        
        Parameters
        ----------
        
        N : int
            Index of state (or site). Signatures make sense only for localized
            states
            
            
        Example
        -------
        
        >>> import quantarhei as qr
        >>> agg = qr.TestAggregate("dimer-2")
        >>> agg.set_coupling_by_dipole_dipole()
        ...
        ... # build and diagonalize
        >>> agg.build()  
        >>> agg.diagonalize()

        >>> agg.get_state_signature_by_index(2)
        ((0, 1), ())
        
        """
        return self.vibsigs[N]
            

def _strip_max_coef(indx, sqrs):
    """Returns the index of the maximum coefficient and the coefficient
    and sets the maximum value to zero.
    
    """
    imax = numpy.argmax(sqrs)
    sqr = sqrs[imax]
    
    # by setting the square to -1 we exclude it from ever becoming a max again
    sqrs[imax] = -1.0
    
    return imax, sqr
        