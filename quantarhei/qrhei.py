# -*- coding: utf-8 -*-
import argparse
import subprocess

from quantarhei import Manager

    
def main():

    parser = argparse.ArgumentParser(
            description='Quantarhei Package Driver')
 
    parser.add_argument("script", metavar='script', type=str,
                        help='script file to be processed')    
    parser.add_argument("-v", "--version", action="store_true",
                        help="shows Quantarhei package version")
    parser.add_argument("-i", "--info", action='store_true', 
                        help="shows detailed information about Quantarhei"+
                        " installation")
    parser.add_argument("-s", "--silent", action='store_true', 
                        help="no output from qrhei script itself")
    parser.add_argument("-p", "--parallel", action='store_true', 
                        help="executes the code in parallel")
    
    args = parser.parse_args() 
    
    flag_parallel = args.parallel
    flag_silent = args.silent

    #
    # show longer info
    #
    if args.info:
        print("")
        print("qrhei: Quantarhei Package Driver")
        print("")
        print("MPI parallelization enabled: ", flag_parallel)
        if not args.version:
            print("Quantarhei package version: ", Manager().version)
            
    #
    # show just Quantarhei version number
    #
    if args.version:
        print("Quantarhei package version: ", Manager().version)

        
    #
    # Running a script
    #
    scr = args.script
    
    if not flag_silent:
        print("Running Quantarhei (python) script file: ", scr)
        print(" --- output below ---")
    
    if not flag_parallel:
        # running the script within the same interpreter
        exec(open(scr).read(), globals())
    else:
        # running MPI with proper parallel configuration
        # FIXME: I want "qrhei" as the running process instead of "python"
        p = subprocess.Popen('mpirun -n 4 python '+scr,
                             shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        for line in iter(p.stdout.readline, ''):
        #for line in p.stdout.readlines():
            ln = line.decode()
            # line is returned with a \n character at the end 
            #ln = ln[0:len(ln)-2]
            print(ln, end="", flush=True)
        retval = p.wait()    
        
        print(retval)
    
