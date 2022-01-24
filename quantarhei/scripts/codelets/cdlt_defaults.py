#
#
###############################################################################
#
# QTask Codelet: dafaults
#
###############################################################################
#
#
# Usage
# -----
# 
# defaults:
#     stdio: True/False
#     files_to_stdio: False/True
#     
#



if INP.defaults:
    
    _default_timeaxis = qr.TimeAxis()
    _default_energy_units = set_or_default(INP.defaults["energy_units"],"1/cm")
    _default_basis = set_or_default(INP.defaults["basis"],"sites")
    
else:

    _default_energy_units = "1/cm"
    _default_basis = "sites"   