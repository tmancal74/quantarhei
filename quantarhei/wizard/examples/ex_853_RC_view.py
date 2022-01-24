# -*- coding: utf-8 -*-
import os
import quantarhei as qr

input_file = "ex_853_RC.yaml"
INP = qr.Input(input_file, show_input=False)

sing = (INP.disorder and (INP.N_realizations == 1))

if INP.single_realization or sing:
    
    print("Single realization 2D map")
    name1 = "cont_p_re"
    
    append_to_dirname = INP.append_to_dirname
    
    file = os.path.join("sim_"+INP.location_of_vibrations
                        +append_to_dirname,name1+".qrp")
    
    print("Loading file:",file)
    
    cont = qr.load_parcel(file)
    sp = cont.get_spectrum(0)
    
    sp.plot(show=False)
    sp.savefig("plot.png")
else:
    
    if INP.disorder:
        
        if INP.N_realizations:
        
            print("Single realization 2D map")
       
        else:
            pass

 
