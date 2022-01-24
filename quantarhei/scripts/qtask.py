# -*- coding: utf-8 -*-

import argparse
import traceback

import quantarhei as qr


def code_hamiltonian(INP, status):
    
    code = "# Ahoj"
    
    return code
    


def tasks(tlist, INP):
    """Go through the list of tasks
    
    """
    
    status = []
    code = ""

    if isinstance(tlist, str):
        tlist = tlist.split()
    print("Identified ", len(tlist), "tasks")
    for task in tlist:
        print("Processing: ", task)
        
        if task == "hamiltonian":
            code += code_hamiltonian(INP, status)
            
    return code
        

def do_process(args):
    """Process arguments of the script
    
    """
    
    if args.filename:
        fname = args.filename[0]
        print("Processing file: ", fname)
    
        try:
            INP = qr.Input(fname)
            code = tasks(INP.tasks, INP)
        except:
            print(traceback.format_exc())
            return
        
    if args.output:
        
        print("Script saved as: ", args.output)
        
    else:
        
        print("\nResulting code:\n")
        print(code)

    

def main():
    
    global parser_list
    global parser_fetch
    
    parser = argparse.ArgumentParser(
            description='Quantarhei Task Compiler')
    
    
    #
    # Driver options
    #
    parser.add_argument("-o", "--output", metavar="output", type=str, 
                        help="output file")
    parser.add_argument("filename", metavar='filename', type=str, 
                        help='task file to be processed', nargs=1)    
    parser.set_defaults(func=do_process)
    
    #
    # Parsing all arguments
    #
    args = parser.parse_args()  

#    #
#    # show longer info
#    #
#    if args.info:
#        qr.printlog("\n" 
#                   +"qtask: Quantarhei Task Compiler\n",
#                   verbose=True, loglevel=0)
##                   +"\n"
##                   +"MPI parallelization enabled: ", flag_parallel,
##                    verbose=True, loglevel=0)
#        if not args.version:
#            qr.printlog("Package version: ", qr.Manager().version, "\n",
#                  verbose=True, loglevel=0)
#        return
#            
#    #
#    # show just Quantarhei version number
#    #
#    if args.version:
#        qr.printlog("Quantarhei package version: ", qr.Manager().version, "\n",
#                  verbose=True, loglevel=0)
#        return
#    
    print("\nQTask: Quantarhei Task Compilator\n")
        
    try:      
        if args.func:
            args.func(args)
    except:
        parser.error("No arguments provided")