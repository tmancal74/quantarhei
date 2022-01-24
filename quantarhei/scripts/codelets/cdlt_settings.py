#
#
###############################################################################
#
# QTask Codelet: settings
#
###############################################################################
#
#
# Usage
# -----
# 
# settings:
#     stdio: True/False
#     files_to_stdio: False/True
#     
#

def set_or_default(a, b):
    """If a is None, default of b is returned instead
    
    """
    if a is not None:
        return a
    else:
        return b
    

if INP.settings:
    
    _stdio = set_or_default(INP.settings["stdio"], True)
    _files_to_stdio = set_of_default(INP.files_to_stdio, False)
    
# defaults
else:
    
    # default is to write reports on standard output
    _stdio = True
    # required files will not be copied on the standard output
    _files_to_stdio = False

