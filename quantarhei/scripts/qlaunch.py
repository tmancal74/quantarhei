# -*- coding: utf-8 -*-
"""
    Quantarhei job launcher
    
    This script is ment to launch Quantarhei jobs on remote machines
    
    The script transfers simulation inputs to the remote machine, launches
    the simulation, monitors it, and transfers the results back to the machine
    from which the job was launched. The simulation inputs are stored in a 
    single directory denoted by a suffix .in. This directory includes input 
    data and configuration settings. It is transferred to the target machine
    and after the simulation is done it is returned back with the suffix .out.


"""

import argparse
import subprocess
import os
import fnmatch
import traceback
import pkg_resources
import sys

import quantarhei as qr



def do_launch():
    
    
    # check for default configuration file: qlaunch.conf
    # check for configuration file within launch directory: DIR/qlaunch.conf
    #    at least one of them has to be provided
    pass


def main():
    
    parser = argparse.ArgumentParser(
            description='Quantarhei Remote Launcher')
    
    parser.add_argument("directory", metavar='directory', type=str, 
                          help='job directory to launch', nargs='?')

    #
    # Driver options
    #
    parser.add_argument("-v", "--version", action="store_true",
                        help="shows Quantarhei package version")
    parser.add_argument("-i", "--info", action='store_true', 
                        help="shows detailed information about Quantarhei"+
                        " installation")
    
    parser.set_defaults(func=do_launch)   
    
    #
    # Parsing all arguments
    #
    args = parser.parse_args() 
    
    
    if len(sys.argv) < 2:
        parser.print_usage()
        qr.exit()
    
    #
    # show longer info
    #
    if args.info:
        qr.printlog("\n" 
                   +"qrhei: Quantarhei Package Driver\n",
                   verbose=True, loglevel=1)
#                   +"\n"
#                   +"MPI parallelization enabled: ", flag_parallel,
#                    verbose=True, loglevel=0)
        if not args.version:
            qr.printlog("Package version: ", qr.Manager().version, "\n",
                  verbose=True, loglevel=1)
        return
            
    #
    # show just Quantarhei version number
    #
    if args.version:
        qr.printlog("Quantarhei package version: ", qr.Manager().version, "\n",
                  verbose=True, loglevel=1)
        return