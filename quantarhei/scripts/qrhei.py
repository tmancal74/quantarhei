# -*- coding: utf-8 -*-
"""
    Driver script for Quantarhei package
    
    
    Author: Tomas Mancal, Charles University, Prague, Czech Republic
    email: mancal@karlov.mff.cuni.cz


"""
import argparse
import subprocess
from pathlib import Path
import os
import sys
import fnmatch

import pkg_resources

import quantarhei as qr


def do_command_run(args):
    """Runs a script 
    
    
    """
    
    m = qr.Manager().log_conf
    
    m.verbosity = args.verbosity

    nprocesses = args.nprocesses
    flag_parallel = args.parallel
    flag_silent = args.silent
    
    if args.silent:
        m.verbosity = 0
        
    #
    # Run benchmark
    #
    if args.benchmark > 0:
        import time

        qr.printlog("Running benchmark no. ", args.benchmark, verbose=True,
                    loglevel=1)
        import quantarhei.benchmarks.bm_001 as bm        
        t1 = time.time()
        bm.main()
        t2 = time.time()
        qr.printlog("... done in", t2-t1, "sec", verbose=True,
                    loglevel=1)
        
        return

    
    #
    # Script name
    # 
    if args.script:
        scr = args.script[0]

    #
    # Greeting 
    #
    qr.printlog("Running Quantarhei (python) script file: ", scr,
                verbose=True, loglevel=3)
        
    
    #
    # Run serial or parallel 
    #
        
    if flag_parallel:
        
        #
        # get parallel configuration
        #
        cpu_count = 0
        try:
            import multiprocessing
            cpu_count = multiprocessing.cpu_count()
        except (ImportError, NotImplementedError):
            pass        
        
        prl_exec = "mpirun"
        prl_n = "-n"
        
        if cpu_count != 0:
            prl_np = cpu_count
        else:
            prl_np = 4
            
        if nprocesses != 0:
            prl_np = nprocesses
        
        engine = "qrhei -s "
        
        # running MPI with proper parallel configuration
        prl_cmd = prl_exec+" "+prl_n+" "+str(prl_np)+" "
        cmd = prl_cmd+engine+scr
        if not flag_silent:
            print("System reports", cpu_count,"processors")
            print("Starting parallel execution with",prl_np,
            "processes (executing command below)")
            print(cmd)
            print("")
        p = subprocess.Popen(cmd,
                             shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)

        if not flag_silent:
            print(" --- output below ---")
        
        # read and print output
        for line in iter(p.stdout.readline, b''):
        #for line in p.stdout.readlines():
            ln = line.decode()
            # line is returned with a \n character at the end 
            # ln = ln[0:len(ln)-2]
            print(ln, end="", flush=True)
            
        retval = p.wait()    
        
    else:
        
        qr.printlog(" --- output below ---", verbose=True, loglevel=0)
        # running the script within the same interpreter
        exec(open(scr).read(), globals())
        
        retval = 0        
        
    #
    # Saying good bye
    #
    if retval == 0:
        qr.printlog("", verbose=True, loglevel=0)
        qr.printlog(" --- output above --- ", verbose=True, loglevel=0)
        qr.printlog("Finished sucessfully; exit code: ", retval,
                    verbose=True, loglevel=0)
    else:
        qr.printlog("Warning, exit code: ", retval, verbose=True, loglevel=0)
        

def do_command_test(args):
    """Runs Quantarhei tests
    
    """
    
    qr.printlog("Running tests", loglevel=0)


def _match_filenames(filenames, pattern, add_stars=False):
    
    if add_stars:
        if not pattern.startswith("*"):
            pattern = "*"+pattern
        if not pattern.endswith("*"):
            pattern = pattern+"*"
            
    return fnmatch.filter(filenames, pattern) 
    
   
def do_command_list(args):
    """Lists files for Quantarhei
    
    """
    global parser_list
    
    if args.examples:
        qr.printlog("Listing available examples ...", loglevel=0)
    
        import quantarhei.wizard.examples as exmpl
         
        filenames = exmpl._available_examples
        
        if args.glob:
            pattern = args.glob
            matching = _match_filenames(filenames, pattern, add_stars=True)
        else:
            matching = filenames
            
        for ex in matching:
            qr.printlog("    "+ex, loglevel=0)

    else:
        parser_list.print_help()


def do_command_fetch(args):
    """Fetches files for Quantarhei
    
    """

    global parser_fetch

    def get_number(file):
        """Returns the number part of the file name
        
        """
        parts = file.split(sep="_")
        return parts[1]
    
    if args.examples:
        qr.printlog("Fetching example(s) ...", loglevel=0)
    
        import quantarhei.wizard.examples as exmpl
         
        filenames = exmpl._available_examples
        data_files = exmpl._available_data
        
        if args.glob:
            pattern = args.glob
            matching = _match_filenames(filenames, pattern, add_stars=True)
        else:
            matching = []
        
        if len(matching) > 0:
            newmatch = matching.copy()
            for match in matching:
                data_pattern = get_number(match)
                newmatch += _match_filenames(data_files, data_pattern,
                                             add_stars=True)
            matching = newmatch 
            
        if len(matching) > 0:           
            for filename in matching:
                
                resource_package = "quantarhei"
                resource_path = '/'.join(('wizard', 'examples', filename))

                content = pkg_resources.resource_string(resource_package,
                                                        resource_path)
                
                over = True

                if os.path.isfile(filename):
                    qr.printlog("File", filename, "already exists!")
                    answr = input(" Overwrite? (y/n) [n]")
                    if answr == "y":
                        over = True
                    else:
                        over = False
                elif os.path.isdir(filename):
                    qr.printlog("Directory with the name", filename,
                                "already exists!")
                    qr.printlog("Aborting fetching")
                    over = False
                    
                if over:
                    with open(filename, "w") as file:
                        file.write(content.decode("utf-8"))
                    
                    qr.printlog("    "+filename, loglevel=0)

        else:
            qr.printlog("No matching examples found", loglevel=0)

    else:
        parser_fetch.print_help()


def do_command_config(args):
    """Configures Quantarhei
    
    """
    
    qr.printlog("Setting configuration", loglevel=0)


def do_command_report(args):
    """Reports on Quantarhei and the system
    
    """
    
    qr.printlog("Probing system configuration", loglevel=0)

    
def main():
    
    global parser_list
    global parser_fetch
    
    parser = argparse.ArgumentParser(
            description='Quantarhei Package Driver')
    
    
    subparsers = parser.add_subparsers(help="Subcommands")

    
    #
    # Driver options
    #
    parser.add_argument("-v", "--version", action="store_true",
                        help="shows Quantarhei package version")
    parser.add_argument("-i", "--info", action='store_true', 
                        help="shows detailed information about Quantarhei"+
                        " installation")
    parser.add_argument("-y", "--verbosity", type=int, default=5, 
                        help="defines verbosity between 0 and 10")
 
    
    #
    # Subparser for command `run`
    #
   
    parser_run = subparsers.add_parser("run", help="Script runner")
    
    parser_run.add_argument("script", metavar='script', type=str, 
                          help='script file to be processed', nargs=1)
    parser_run.add_argument("-s", "--silent", action='store_true', 
                          help="no output from qrhei script itself")
    parser_run.add_argument("-p", "--parallel", action='store_true', 
                          help="executes the code in parallel")
    parser_run.add_argument("-n", "--nprocesses", type=int, default=0,
                          help="number of processes to start")
    parser_run.add_argument("-b", "--benchmark", type=int, default=0, 
                          help="run one of the predefined benchmark"
                          +"calculations")
    
    parser_run.set_defaults(func=do_command_run)
    
    #
    # Subparser for command `test`
    #

    parser_test = subparsers.add_parser("test", help="Test runner")
    
    parser_test.set_defaults(func=do_command_test)    
    
    #
    # Subparser for command `fetch`
    #

    parser_fetch = subparsers.add_parser("fetch", help="Fetches examples,"
                                        +" benchmarks, tutorials, templates"
                                        +" and configuration files")

    parser_fetch.add_argument("glob", metavar='glob', type=str, 
                              help='file name', nargs="?")
    parser_fetch.add_argument("-e", "--examples", action='store_true',
                              help="fetches a specified example file")
    
    parser_fetch.set_defaults(func=do_command_fetch)    

    #
    # Subparser for command `list`
    #

    parser_list = subparsers.add_parser("list", help="Lists examples,"
                                    +" benchmarks, tutorials and templates")

    parser_list.add_argument("glob", metavar='glob', type=str, 
                             help='file name', nargs="?")
    parser_list.add_argument("-e", "--examples", action='store_true', 
                             help="list all available example files")
    
    parser_list.set_defaults(func=do_command_list)    


    
    #
    # Subparser for command `config`
    #

    parser_conf = subparsers.add_parser("config", help="Configures Quantarhei")
    
    parser_conf.set_defaults(func=do_command_config)    

    #
    # Subparser for command `report`
    #

    parser_report = subparsers.add_parser("report", help=
                                "Probes Quantarhei as system configurations")
    
    parser_report.set_defaults(func=do_command_report)    
    
    #
    # Parsing all arguments
    #
    args = parser.parse_args()       

    #
    # show longer info
    #
    if args.info:
        qr.printlog("\n" 
                   +"qrhei: Quantarhei Package Driver\n",
                   verbose=True, loglevel=0)
#                   +"\n"
#                   +"MPI parallelization enabled: ", flag_parallel,
#                    verbose=True, loglevel=0)
        if not args.version:
            qr.printlog("Package version: ", qr.Manager().version, "\n",
                  verbose=True, loglevel=0)
        return
            
    #
    # show just Quantarhei version number
    #
    if args.version:
        qr.printlog("Quantarhei package version: ", qr.Manager().version, "\n",
                  verbose=True, loglevel=0)
        return
    
        
    try:      
        if args.func:
            args.func(args)
    except:
        parser.error("No arguments provided")

        
    
