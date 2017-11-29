# -*- coding: utf-8 -*-
import numpy

from .aggregate_spectroscopy import AggregateSpectroscopy
import quantarhei as qr

class AggregateExcitonAnalysis(AggregateSpectroscopy):
    
    
    def get_expansion_squares(self, state=0):
        """Returns the squares of expansion coefficients of an excitonic state
        
        """
        coefs = self.SS[:,state]
        sqrs  = self.SS[:,state]**2
        indx = [i for i in range(self.HH.shape[0])]
        return indx, coefs, sqrs
        
    
    # FIXME: what happens when there are no vibrational modes
    def get_state_signature_by_index(self, N):
        """Return aggregate vibronic state signature  by its index
        
        """
        return self.vibsigs[N]
    
    
    def report_on_expansion(self, state=0, N=5):
        """Prints a short report on the composition of an exciton state
        
        """
        try:
            from terminaltables import AsciiTable
        except:
            raise Exception("Get terminaltables package "
                            +"for this functionality")
        
        #idxs = numpy.zeros(N, dtype=numpy.int)
        #cofs = numpy.zeros(N, dtype=qr.REAL)
    
        # to be a method
        indx, coefs, sqrs = self.get_expansion_squares(state)
    
        table_data = []
        table_data.append(["index","squares", "coefficients",
                           "state signatures"])
        for i in range(N):
            imax, sqr = _strip_max_coef(indx, sqrs)
            coef = coefs[imax]
            # to be a method
            sta = self.get_state_signature_by_index(imax)
            table_data.append([imax, sqr, coef, sta])
            
            #print(imax, "\t", coef,"\t", sta)
            
        table = AsciiTable(table_data)
        print(table.table)


    def get_intersite_mixing(self, state1=0, state2=0):
        """Returns inter site mixing ration
        
        Inter site mixing ratio gives the probability that
        states ``state1`` and ``state2`` are localized
        on different pigments.
        
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

def _strip_max_coef(indx, sqrs):
    """Returns the index of the maximum coefficient and the coefficient
    and sets the maximum value to zero.
    
    """
    imax = numpy.argmax(sqrs)
    sqr = sqrs[imax]
    
    sqrs[imax] = 0.0
    
    return imax, sqr
        