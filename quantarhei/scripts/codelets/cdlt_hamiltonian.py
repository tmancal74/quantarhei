#
#
###############################################################################
#
# QTask Codelet: hamiltonian
#
###############################################################################
#
#
# Usage
# -----
# 
# hamiltonian:
#     matrix:
#         12000.0
#         0.0           12300.0
#         0.0           300.0      12500.0
#     units: 1/cm
#     multiplicity: 2
#     output:
#         basis: eigenstates
#         file: ham.dat

class QTaskException(Exception):
    pass


if INP.hamiltonian:
    
    #
    # Get multiplicity
    #
    if INP.hamiltonian["multiplicity"]:
        mult = INP.hamiltonian["multiplicity"]
    else:
        mult = 1

    #
    # Get units
    #      
    if INP.hamiltonian["units"]:
        units = INP.hamiltonian["units"]
    else:
        units = "int"
        
    
    #
    # Either matrix or file keywords have to be specified
    #
    error = False
    have_data = False
    try:
        data = INP.hamiltonian["matrix"]
        have_data = True
    except KeyError:
        error = True
    
    if error:
        try:
            ham_file = INP.hamiltonian["file"]
            error = False
        except KeyError:
            error = True
    
    if error:
        raise QTaskException("`matrix` or `file` keywords have to be provided")
    
    if have_data:
        H = qr.Hamiltonian(data=data)
    
else:
    
    raise QTaskException("Hamiltonian is not provided")