# -*- coding: utf-8 -*-
def rename_function(ss, oldname, newname):
    """Replaces all occurences of a name by a new name
    
    """
    #return ss.replace("conjugate","numpy.conj")
    return ss.replace(oldname, newname)
    
def fce2array(sr, pat):
    """Converts functions into arrays 
    
    """
    se = "".join(sr)
    so = ""
    ln = len(se)
    while ln > 0:
        # find pattern
        pos = se.find(pat)
        if pos < 0:
            break
        # position just behind the pattern
        pos += len(pat)
        sl = list(se)
        # exchange ( for [
        if sl[pos] == "(":
            sl[pos] = "["
        se = "".join(sl)
        # save everything in front of the pattern
        so += se[0:pos]
        se = se[pos:ln]
        # find clossing braket
        pos2 = se.find(")")
        # echange ) for ]
        sl = list(se)
        if sl[pos2] == ")":
            sl[pos2] = "]"
        se = "".join(sl)
    
        ln = len(se)
    so += se
    return so

def python_code(ss, arrays=None):
    """Generate Python code with numpy functions
    
    """
    sr = rename_function(ss,"conjugate","numpy.conj")    
    sr = rename_function(sr,"exp","numpy.exp")
    if arrays is not None:
        for ar in arrays:    
            sr = fce2array(sr,ar)
    return sr

def fortran_code(ss, arrays=None):
    """Generate Fortran code with numpy functions
    
    """
    sr = rename_function(ss,"conjugate","conjg")    
    #sr = rename_function(sr,"exp","numpy.exp")
    #for ar in arrays:    
    #    sr = fce2array(sr,ar)
    return sr
