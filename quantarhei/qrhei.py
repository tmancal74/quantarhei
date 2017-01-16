# -*- coding: utf-8 -*-
import argparse
import sys
import subprocess

from quantarhei import Manager

def runProcess(exe):    
    p = subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while(True):
      retcode = p.poll() #returns None while subprocess is running
      line = p.stdout.readline()
      yield line
      if(retcode is not None):
        break
    
def process_output(line):
    print(line)
    
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
    
    args = parser.parse_args() 
    
    flag_parallel = False
    flag_silent = args.silent

    # show longer info
    if args.info:
        print("")
        print("qrhei: Quantarhei Package Driver")
        print("")
        print("MPI parallelization enabled: ", flag_parallel)
        if not args.version:
            print("Quantarhei package version: ", Manager().version)
            
    
    # show just Quantarhei version number
    if args.version:
        print("Quantarhei package version: ", Manager().version)

        
    
    scr = args.script
    if not flag_silent:
        print("Running Quantarhei (python) script file: ", scr)
        print(" --- output below ---")
    
    if not flag_parallel:
        # running the script within the same interpreter
        exec(open(scr).read(), globals())
    else:
        # running MPI with proper parallel configuration
        pass
    
    
