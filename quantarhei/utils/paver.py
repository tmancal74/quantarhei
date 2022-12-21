# -*- coding: utf-8 -*-

import subprocess
import sys

def execute_paver(args, util="paver"):
    """Executes 'paver' utility
    
    """
    pargs = ['paver']+args
    process = subprocess.Popen(pargs,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT,
                               universal_newlines=True)

    for line in process.stdout:
        sys.stdout.write(line)
        sys.stdout.flush()

    #while True:
    #    output = process.stdout.readline()
    #    print(output.strip())
    #    # Do something else
    #    return_code = process.poll()
    #    if return_code is not None:
    #        #print(util+' finished with return code', return_code)
    #        # Process has finished, read rest of the output
    #        for output in process.stdout.readlines():
    #            print(output.strip())
    #        break