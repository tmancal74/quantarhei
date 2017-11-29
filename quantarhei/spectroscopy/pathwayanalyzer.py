# -*- coding: utf-8 -*-

from ..core.managers import UnitsManaged
from ..core.wrappers import deprecated

class LiouvillePathwayAnalyzer(UnitsManaged):
    """Class providing methods for Liouville pathway analysis
    
    
    Parameters
    ----------
    
    pathways : list
        List of pathways to analyze
    
    """

    def __init__(self, pathways=None):
        self.pathways = pathways


    def set_pathways(self, pathways):
        """Sets pathways to be analyzed
        
        """
        self.pathways = pathways


    def get_pathways(self):
        """Returns Liouville pathways
        
        """
        return self.pathways
    
    @deprecated
    def max_pref(self, pathways):
        """Return the maximum of pathway prefactors
        
        
        Parameters
        ----------
        
        pathways : list
            List of Liouville pathways
            
            
        Returns
        -------
        
        pmax : float
            Maximum prefactor of the pathways
            
        rec : int
            position of the pathway with the maximum prefactor
        
        """
        
        if pathways is None:
            pthways = self.pathways
        else:
            pthways = pathways
        
        pmax = 0.0
        k = 0
        rec = -1
        for pway in pthways:
            if pway.pref > pmax:
                rec = k
                pmax = pway.pref
            k += 1
            
        return (pmax, rec)

    def max_amplitude(self):
        """Return the maximum of pathway prefactors
        
        
        Parameters
        ----------
        
        pathways : list
            List of Liouville pathways
            
            
        Returns
        -------
        
        pmax : float
            Maximum prefactor of the pathways
            
        rec : int
            position of the pathway with the maximum prefactor
        
        """
        return max_amplitude(self.pathways)

    def select_pref_GT(self, val, pathways=None, replace=True,
                       verbose=False):
        """Select all pathways with prefactors greater than a value
        
        """
        
        if pathways is None:
            pthways = self.pathways
        else:
            pthways = pathways
            
        selected = []
        for pway in pthways:
            if numpy.abs(pway.pref) > val:
                selected.append(pway)
        
        if verbose:
            print("Selected", len(selected), "pathways")
            
        # if the pathways were not specified from argument then return selected
        # to self
        if (pathways is None) and replace:
            self.pathways = selected
            
        return selected

    
    def select_frequency_window(self, window, pathways=None, replace=True, 
                                verbatime=False):
        """Selects pathways with omega_1 and omega_3 in a certain range
        
        """
        if pathways is None:
            pthways = self.pathways
        else:
            pthways = pathways

        om1_low = self.convert_energy_2_internal_u(window[0])
        om1_upp = self.convert_energy_2_internal_u(window[1])
        om3_low = self.convert_energy_2_internal_u(window[2])
        om3_upp = self.convert_energy_2_internal_u(window[3])
        
        
        selected = []
        
        for pway in pthways:
            ne = len(pway.frequency)
            #om1 = numpy.abs(pway.frequency[0])
            om1 = numpy.abs(pway.get_interval_frequency(0))
            #om3 = pway.frequency[ne-2]
            om3 = numpy.abs(pway.get_interval_frequency(ne-2))
            
            if (((om1 >= om1_low) and (om1 <= om1_upp)) and 
                ((om3 >= om3_low) and (om3 <= om3_upp))):
                selected.append(pway) 
                
        if verbatime:
            print("Selected", len(selected), "pathways")

        if (pathways is None) and replace:
            self.pathways = selected
            
        return selected


    def select_omega2(self, interval, pathways=None, replace=True,
                      verbatime=False):
        """Selects pathways with omega_2 in a certain interval
        
        """    
        if pathways is None:
            pthways = self.pathways
        else:
            pthways = pathways

        om2_low = self.convert_energy_2_internal_u(interval[0])
        om2_upp = self.convert_energy_2_internal_u(interval[1])  
        
        selected = []
        
        for pway in pthways:
            ne = len(pway.frequency)
            #om2 = pway.frequency[ne-3]
            om2 = pway.get_interval_frequency(ne-3)
            
            if (om2 >= om2_low) and (om2 <= om2_upp):
                selected.append(pway)
                
        if verbatime:
            print("Selected", len(selected), "pathways")

        if (pathways is None) and replace:
            self.pathways = selected
        
        return selected
    
    
    def order_by_pref(self, pthways):
        """Orders the list of pathways by pathway prefactors
        
        """
        lst = sorted(pthways, key=lambda pway: abs(pway.pref), reverse=True)
        return lst
    

def max_amplitude(pathways):
    """Return the maximum of pathway prefactors
    
    
    Parameters
    ----------
    
    pathways : list
        List of Liouville pathways
        
        
    Returns
    -------
    
    pmax : float
        Maximum prefactor of the pathways
        
    rec : int
        position of the pathway with the maximum prefactor
    
    """
    
    pmax = 0.0
    k = 0
    rec = -1
    for pway in pathways:
        if pway.pref > pmax:
            rec = k
            pmax = pway.pref
        k += 1
        
    return (pmax, rec)
