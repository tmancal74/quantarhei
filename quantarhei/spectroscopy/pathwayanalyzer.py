# -*- coding: utf-8 -*-
"""Liouville pathway analysis module

   This module defines classes and function to help in 
   Liouvile pathway analysis
   
   Classes
   -------
   
   LiouvillePathwayAnalyzer
      
   
   Functions
   ---------
   
   max_amplitude(pathways)
   
   
   
   

"""
import os
import numpy

from ..core.managers import UnitsManaged, Manager
from ..core.wrappers import deprecated

from ..core.units import cm2int
from ..core.units import convert
from ..core.parcel import load_parcel
from .. import REAL, COMPLEX

from ..core.time import TimeAxis
from ..core.dfunction import DFunction



class LiouvillePathwayAnalyzer(UnitsManaged):
    """Class providing methods for Liouville pathway analysis
    
    Parameters
    ----------
    
    pathways : list
        List of pathways to analyze
    
    
    Methods
    -------
    
    set_pathways(pathways=None)
        Set the pathways for the analysis
        
    get_pathways()
        Returns current list of pathways
        
    max_amplitude()
        Returns a tuple containing the value of the amplitude of the pathway
        from current list of pathways with the maximum amplitude and its index
        in the list
        
    select_amplitide_GT(val, replace=True, verbose=False)
        In the current list of pathways this method selects those that have 
        absolute value of amplitude larger that a given value
        
    
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
    
    def get_number_of_pathways(self):
        """Returns the number of available pathways
        
        """
        return(len(self.pathways))
    
    
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


    @deprecated
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


    def select_amplitude_GT(self, val, replace=True,
                       verbose=False):
        """Select all pathways with abs value of prefactors greater than a value
    
        """
        selected =  select_amplitude_GT(val, self.pathways, verbose=verbose)
        # if the pathways were not specified from argument then return selected
        # to self
        if replace:
            self.pathways = selected
        else: 
            return selected
    

    def select_frequency_window(self, window, replace=True, 
                                verbose=False):
        """Selects pathways with omega_1 and omega_3 in a certain range
        
        """
        selected = select_frequency_window(window, self.pathways, verbose)

        if replace:
            self.pathways = selected
        else:
            return selected


    def select_omega2(self, interval, replace=True, verbose=False):
        """Selects pathways with omega_2 in a certain interval
        
        """    
        selected = select_omega2(interval, self.pathways, verbose)

        if replace:
            self.pathways = selected
        else:
            return selected
    
    @deprecated
    def order_by_pref(self, pthways):
        """Orders the list of pathways by pathway prefactors
        
        """
        lst = sorted(pthways, key=lambda pway: abs(pway.pref), reverse=True)
        return lst


    def order_by_amplitude(self, replace=True):
        
        orderred = order_by_amplitude(self.pathways)
        if replace:
            self.pathways = orderred
        else:
            return orderred


    def select_sign(self, sign, replace=True):
        
        selected = select_sign(self.pathways, sign)
        if replace:
            self.pathways = selected
        else:
            return selected
        
    def select_type(self, ptype="REPH", replace=True):
        
        selected = select_type(self.pathways, ptype)
        if replace:
            self.pathways = selected
        else:
            return selected
    
        return selected

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


def select_amplitude_GT(val, pathways, verbose=False):
    """Select all pathways with abs value of prefactors greater than a value
    
    Parameters
    ----------
    
    val : float
        Value to compare with the pathway prefactor
        
    pathways : list of pathways
        List of the Liouville pathways to analyze
        
    verbose: bool
        If set True, number of selected pathways is reported
        
        
    Returns
    -------
    
    selected : list
        List of selected pathways
        
    """
    
    pthways = pathways
        
    selected = []
    for pway in pthways:
        if numpy.abs(pway.pref) > val:
            selected.append(pway)
    
    if verbose:
        print("Selected", len(selected), "pathways")
        
        
    return selected


def select_frequency_window(window, pathways, verbose=False):
    """Selects pathways with omega_1 and omega_3 in a certain range
        
    """

    pthways = pathways
    m = Manager()

    om1_low = m.convert_energy_2_internal_u(window[0])
    om1_upp = m.convert_energy_2_internal_u(window[1])
    om3_low = m.convert_energy_2_internal_u(window[2])
    om3_upp = m.convert_energy_2_internal_u(window[3])
    
    
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
            
    if verbose:
        print("Selected", len(selected), "pathways")
        
    return selected
 

def select_omega2(interval, pathways, secular=True, 
                  tolerance=10.0*cm2int, verbose=False):
    """Selects pathways with omega_2 in a certain interval
    
    """    

    pthways = pathways
    m = Manager()

    om2_low = m.convert_energy_2_internal_u(interval[0])
    om2_upp = m.convert_energy_2_internal_u(interval[1])  
    
    selected = []
    
    for pway in pthways:
        ne = len(pway.frequency)
        #print("Number of frequencies: ", ne)
        om2 = pway.get_interval_frequency(ne-3)
        
        # pathways with transfer
        if ne > 4:
            # check previous frequency (feeding)
            om2_2 =  pway.get_interval_frequency(ne-4)
        
            if secular:
                # case with feeding frequency small
                if numpy.abs(om2_2) <= tolerance:
                    if (om2 >= om2_low) and (om2 <= om2_upp):
                        selected.append(pway) 
                        #print("selected")
                
                # case of fast feeding frequency
                elif numpy.abs(om2 - om2_2) <= tolerance:
                    if (om2 >= om2_low) and (om2 <= om2_upp):
                        selected.append(pway)
                        #print("selected")
                        
            else:
               if (om2 >= om2_low) and (om2 <= om2_upp):         
                   selected.append(pway)
                   #print("selected")
          
        else:
            if (om2 >= om2_low) and (om2 <= om2_upp):
                selected.append(pway)
                #print("selected")
                
            
    if verbose:
        print("Selected", len(selected), "pathways")
        
    return selected


#def select_incoherent(pathways):
#    
#    for pway in pathways:
        

def order_by_amplitude(pthways):
    """Orders the list of pathways by pathway prefactors
    
    """
    lst = sorted(pthways, key=lambda pway: abs(pway.pref), reverse=True)
    return lst


def select_sign(pathways, sign):
    """Selects all pathways depending on the overall sign """
    selected = []
    
    pos = False
    if sign > 0.0:
        pos = True
        
    for pway in pathways:
        if pos:
            if pway.sign > 0.0:
                selected.append(pway)
        else:
            if pway.sign < 0.0:
                selected.append(pway)
    
    return selected


def select_type(pathways, stype):
    """Selects all pathways of a given type
    
    Parameters
    ----------
    
    pathways : list
        List of pathways to be analyzed
        
    stype : str
        Type of the pathways to be selected. 
        Possible values are "REPH", "NONR".
        
    """
    di = dict(REPH="R", NONR="NR")
    selected = []
    
    for pw in pathways:
        if pw.pathway_type == di[stype]:
            selected.append(pw)
    
    return selected
    

def select_by_states(pathways, states):
    """Returns one pathway which goes through a given pattern of states 
    
    Returns unique pathway which goes through a given pattern of states
    or does not return anything

    The following diagram 

            |12       12|  
        --->|-----------|  
            |           |  
            |48       12|      2.0
        --->|-----------|  
            |           |  
            |21       12|      0.0
          >>|***********|<< 
            |           | 
            |6         6|      0.0
        --->|-----------|  
            |           |  
            |0         6|      -2.0
            |-----------|<---  
            |           |  
            |0         0|  
    
    turns into the following tuple
    
    ((0,0), (0,6), (6,6), (21, 12), (48, 12), (12,12))
    
    
    Parameters
    ----------
    
    pathways : list, tuple
        List of tuple of states to analyze
        
    states : list, tuple
        List of dyads which describe the states involved in a given Liouville
        pathway. States are listed from left to right, from bottom to the top
        
    
    
    """

    # loop over all pathways
    for pw in pathways:

        ch = -1
        # loop over the Feynman diagram
        for k in range(pw.states.shape[0]):
            
            # check if prescribed states match the pathway states
            if not ((pw.states[k, 0] == states[k+1][0]) &
                    (pw.states[k, 1] == states[k+1][1])):
                # leave if there is a mismatch
                break
            
            ch += 1
            
        # if all lines match, return the pathway
        if ch == k:
            return pw
        

def look_for_pathways(name="pathways", ext="qrp",
                      check=False, directory="."):
    """Load pathways by t2
    
    """
    
    # get all files with the name_*.ext pattern
    import glob
    import os.path
    
    path = os.path.join(directory,name+"_*."+ext)
    files = glob.glob(path)

    # get a list of t2s
    t2s = []
    for fl in files:
        t2 = float(fl.split("_")[1].split("."+ext)[0])
        t2s.append(t2)
    
    t2s = numpy.array(t2s)
    t2s = numpy.sort(t2s)
    # if check == True load them all and check they are list of pathways
    if check:
        pass
    
    
    return t2s


def load_pathways_by_t2(t2, name="pathways", ext="qrp", directory=".",
                        tag_type=REAL):
    """Load pathways by t2 time
    
    """
    
    t2_str = str(t2)
    fname = name+"_"+t2_str+"."+ext
    path = os.path.join(directory, fname)
    try:
        pw = load_parcel(path)
    except:
        print("Error while loading")
        return []
    
    return pw


def save_pathways_by_t2(t2, name="pathways", ext="qrp", directory=".",
                        tag_type=REAL):
    pass


# FIXME: This is done in an inefficient way, loading the files multiply
def get_evolution_from_saved_pathways(states, name="pathways", ext="qrp", 
                                      directory=".", tag_type=REAL, repl=0.0):
    """Reconstructs the evolution of the pathway contribution in t2 time
    
    """
    t2s = look_for_pathways(name=name, ext=ext, directory=directory)
   

    # order t2s
    t2s = numpy.sort(t2s)
    evol = _get_evol(t2s, states, name, ext, directory, repl=repl)
    
    dt = t2s[1] - t2s[0]
    length = len(t2s)
    
    taxis = TimeAxis(t2s[0], length, dt)
    
    ii = 0
    for tt in t2s:
        if tt != taxis.data[ii]:
            raise Exception("The set of available times"+
                            " does not correspond to a continuous time axis")
        ii += 1
    
    return DFunction(x=taxis, y=evol)
    
    #return t2s, evol


def get_prefactors_from_saved_pathways(states, name="pathways", ext="qrp", 
                                      directory=".", tag_type=REAL, repl=0.0):
    """Reconstructs the evolution of the pathway contribution in t2 time
    
    """
    t2s = look_for_pathways(name=name, ext=ext, directory=directory)
   

    # order t2s
    t2s = numpy.sort(t2s)
    evol = _get_pref(t2s, states, name, ext, directory, repl=repl)
    
    dt = t2s[1] - t2s[0]
    length = len(t2s)
    
    taxis = TimeAxis(t2s[0], length, dt)
    
    ii = 0
    for tt in t2s:
        if tt != taxis.data[ii]:
            raise Exception("The set of available times"+
                            " does not correspond to a continuous time axis")
        ii += 1
    
    return DFunction(x=taxis, y=evol)
         

def get_TwoDSpectrum_from_saved_pathways(t2, t1axis, t3axis, name="pathways",
                                         ext="qrp", directory=".",
                                         tag_type=REAL):
    """Returns a 2D spectrum calculated based on the saved Liouville pathways
    
    """
    pwt2 = load_pathways_by_t2(t2, name=name, ext=ext,
                               directory=directory, tag_type=tag_type)
    
    # calculate 2D spectrum from the loaded pathways
    twod = get_TwoDSpectrum_from_pathways(pwt2, t1axis, t3axis)
    
    # return it
    return twod


def get_TwoDSpectrum_from_pathways(pathways, t1axis, t3axis):
    """Returns a 2D spectrum calculated based on submitted Liouville pathways
    
    """    
    from .mocktwodcalculator import MockTwoDSpectrumCalculator
    t2axis = TimeAxis(0.0,1,1.0)
    mcalc = MockTwoDSpectrumCalculator(t1axis, t2axis, t3axis)
    mcalc.bootstrap(rwa = convert(12000.0,"1/cm","int"), pathways=pathways)
    twod = mcalc.calculate()

    return twod


def get_TwoDSpectrumContainer_from_saved_pathways(t1axis, t3axis,
                                                  name="pathways", ext="qrp",
                                                  directory=".",
                                                  tag_type=REAL):
    """Returns a container with 2D spectra calculated from saved pathways
    
    """
    from .twodcontainer import TwoDSpectrumContainer
    
    t2s = look_for_pathways(name=name, ext=ext, directory=directory)
    

    # order t2s
    t2s = numpy.sort(t2s)
    time2 = TimeAxis(t2s[0], len(t2s), t2s[1]-t2s[0])
    tcont = TwoDSpectrumContainer(t2axis=time2)
    tcont.use_indexing_type(itype=time2)
    
    for t2 in t2s:
        pwt2 = load_pathways_by_t2(t2, name=name, ext=ext,
                                   directory=directory, tag_type=tag_type)
        
        
        # calculate 2D spectrum from saved pathways
        twod = get_TwoDSpectrum_from_pathways(pwt2, t1axis, t3axis)
        
        # put it into container
        tcont.set_spectrum(twod, tag=t2)
        
    # return the container
    return tcont
    

        



    
def _is_tuple_of_dyads(states):
    """Check if the object is a tuple or list of dyads
    
    """
    
    def _is_a_dyad(dd):
        return len(dd) == 2
    
    ret = False
    for st in states:
        if _is_a_dyad(st):
            ret = True
        else:
            ret = False
            break
        
    return ret


def _get_evol(t2s, states, name, ext, directory, repl=0.0):
    """Return evolution of a single pathway
    
    """
    
    N = 1
    if _is_tuple_of_dyads(states):
        evol = numpy.zeros(len(t2s), dtype=COMPLEX)
    else:
        N = len(states)
        evol = numpy.zeros((len(t2s),N), dtype=COMPLEX)
        
    k = 0
    for t2 in t2s:
        pws = load_pathways_by_t2(t2, name=name, ext=ext, directory=directory)
        
        if N == 1:
            pw = select_by_states(pws, states)
            if pw is None:
                evol[k] = repl
            else:
                evol[k] = pw.evolfac
        else:
            l = 0
            for st in states:
                pw = select_by_states(pws, st)
                if pw is None:
                    evol[k,l] = repl
                else:
                    evol[k,l] = pw.evolfac
                l += 1
                
        k += 1
    
    return evol

def _get_pref(t2s, states, name, ext, directory, repl=0.0):
    """Return evolution of a single pathway
    
    """
    
    N = 1
    if _is_tuple_of_dyads(states):
        evol = numpy.zeros(len(t2s), dtype=COMPLEX)
    else:
        N = len(states)
        evol = numpy.zeros((len(t2s),N), dtype=COMPLEX)
        
    k = 0
    for t2 in t2s:
        pws = load_pathways_by_t2(t2, name=name, ext=ext, directory=directory)
        
        if N == 1:
            pw = select_by_states(pws, states)
            if pw is None:
                evol[k] = repl
            else:
                evol[k] = pw.pref
        else:
            l = 0
            for st in states:
                pw = select_by_states(pws, st)
                if pw is None:
                    evol[k,l] = repl
                else:
                    evol[k,l] = pw.pref
                l += 1
                
        k += 1
    
    return evol
    
    