# -*- coding: utf-8 -*-
"""
    Driver script for Quantarhei package
    
    
    Author: Tomas Mancal, Charles University, Prague, Czech Republic
    email: mancal@karlov.mff.cuni.cz


"""
import argparse
import subprocess
import os
import fnmatch
import traceback
import pkg_resources

import quantarhei as qr


def do_command_run(args):
    """Runs a script 
    
    
    """
    
    m = qr.Manager().log_conf
    dc = qr.Manager().get_DistributedConfiguration()
    
    #
    # analyzing --verbosity option
    #
    verb = args.verbosity
    try:
        vrbint = int(verb)
        m.verbosity = vrbint
        m.fverbosity = vrbint + 2
    except:
        try:
            vrbint = verb.split(",")
            m.verbosity = int(vrbint[0])
            m.fverbosity = int(vrbint[1])
        except:
            raise Exception("Integer or two comma separated integers required"
                            +" for -y/--verbosity option")

    
    # we set the verbosity lower for other than the leading process
    if dc.rank != 0:
        m.verbosity -= 2
        m.fverbosity -= 2
    
    #
    # log into file
    #
    m.log_to_file = args.logtofile
    
    if m.log_to_file:
        m.log_file_name = args.filename
    
    #
    # parallel options
    #
    nprocesses = args.nprocesses
    flag_parallel = args.parallel
    hostfile=""
    if len(args.hostfile) > 0:
        hostfile = args.hostfile

    
    #
    # some other logging option
    #
    flag_silent = args.silent
    flag_quiet = args.quiet
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
    # Input file name  (-i option)
    #
    if args.inputfile:
        input_file = args.inputfile


    #
    # if the file is yaml, look into it to find the script file name
    #
    # get the file extension
    extsplt = os.path.splitext(scr)
    ext = extsplt[1]
    
    #
    # Reading configuration file
    #
    
    # yaml
    if ext in [".yaml", ".yml"]:
        INP = qr.Input(scr)
        
        # this overrides the -i option 
        input_file = scr
        
        # see if the script name is specified, if not use the conf name + .py
        try:
            script = INP.script
        except:
            script = extsplt[0]+".py"
        # see if path is defined, if not use local directory
        try:
            spath = INP.path
        except:
            spath = "."
        
        scr = os.path.join(spath, script)

    #
    # Greeting 
    #
    if not flag_quiet:
        qr.printlog("",verbose=True, loglevel=qr.LOG_URGENT)
        qr.printlog("Running Quantarhei (python) script file: ", scr,
                    verbose=True, loglevel=qr.LOG_URGENT)
    
    #
    # Run serial or parallel 
    #   

        
    if flag_parallel:
        
        #
        # If this is set to True, we use one more processes than processor
        # number to steer other processes
        #
        # Also set the corresponding flag in the parallel module
        #
        use_steerer = False
        if use_steerer:
            nsteerer = 1
        else:
            nsteerer = 0
            
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
        prl_h = ""
        if len(hostfile) > 0:
            prl_h += " --hostfile "+hostfile+" "
        
        if cpu_count != 0:
            prl_np = cpu_count + nsteerer
        else:
            prl_np = 4
            
        if nprocesses != 0:
            prl_np = nprocesses + nsteerer
        
        #
        engine = "qrhei"
        if m.log_to_file:
            engine += " -lf "+m.log_file_name
        engine += " -y "+str(m.verbosity)+","+str(m.fverbosity)+" run -q"
        engine += " -i "+input_file+" "
        
        # running MPI with proper parallel configuration
        prl_cmd = prl_exec+" "+prl_n+" "+str(prl_np)+" "+prl_h
        cmd = prl_cmd+engine+scr
        
        if not flag_silent:
            qr.printlog("System reports", cpu_count,"processors")
            qr.printlog("Starting parallel execution with",prl_np,
            "processes (executing command below)")
            qr.printlog(">>>", cmd)
            qr.printlog("")
            
        try:
            p = subprocess.Popen(cmd,
                             shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)

            if not flag_silent and (not flag_quiet):
                qr.printlog(" --- output below ---\n", verbose=True,
                            loglevel=1)
                
            # read and print output
            for line in iter(p.stdout.readline, b''):
            #for line in p.stdout.readlines():
                ln = line.decode()
                # line is returned with a \n character at the end 
                # ln = ln[0:len(ln)-2]
                print(ln, end="", flush=True)
                
            retval = p.wait()  
            
        except SystemExit:
            
            qr.printlog("", verbose=True, loglevel=qr.LOG_DETAIL)
            qr.printlog(" --- Exited by SystemExit --- ", verbose=True,
                        loglevel=qr.LOG_DETAIL)
            pass
            
    else:
        
        if not flag_silent and (not flag_quiet):
            qr.printlog(" --- output below ---\n", verbose=True, 
                        loglevel=qr.LOG_URGENT)
        
        # running the script within the same interpreter
        try:

            # launch this properly, so that it gives information
            # on the origin of exceptions
            with open(scr,'U') as fp:
                code = fp.read()
                glbs = globals()
                glbs.update(dict(_input_file_=input_file))
                    
            exec(compile(code, scr, "exec"), glbs)
            
            
        except SystemExit:
        
            qr.printlog("", verbose=True, loglevel=qr.LOG_DETAIL)
            qr.printlog(" --- Exited by SystemExit --- ", verbose=True,
                        loglevel=qr.LOG_DETAIL) 
            
        except: 
            
            print(traceback.format_exc())
        
        retval = 0        
        
    #
    # Saying good bye
    #
    if retval == 0:
        if not flag_silent and (not flag_quiet):
            qr.printlog("", verbose=True, loglevel=1)
            qr.printlog(" --- output above --- \n", verbose=True, 
                        loglevel=qr.LOG_URGENT)
            qr.printlog("Finished sucessfully; exit code: ", retval,
                        verbose=True, loglevel=1)
    else:
        qr.printlog("Warning, exit code: ", retval, verbose=True, 
                    loglevel=qr.LOG_URGENT)
        
    
    

def do_command_test(args):
    """Runs Quantarhei tests
    
    """
    
    qr.printlog("--- Running tests ---", loglevel=qr.LOG_URGENT)
    qr.printlog("")
    qr.printlog("Testing parallel capabilities:", loglevel=qr.LOG_URGENT)
    
    #
    # number of processors
    #
    cpu_count = 0
    try:
        import multiprocessing
        cpu_count = multiprocessing.cpu_count()
    except (ImportError, NotImplementedError):
        pass
    qr.printlog("Processor count:", cpu_count)

    #
    # mpi4py
    #
    have_mpi = False
    try:
        from mpi4py import MPI
        have_mpi = True
    except:
        pass
    qr.printlog("mpi4py installed:", have_mpi)

    #
    # MPI implementation
    #
    cmd = "mpirun -h"
    mpirun = _try_cmd(cmd)
        
    qr.printlog("mpirun available:", mpirun)

    #
    # MPI implementation
    #
    cmd = "mpiexec -h"
    mpiexec = _try_cmd(cmd)
        
    qr.printlog("mpiexec available:", mpiexec)


def _try_cmd(cmd):
    ret = False
    try:
        p = subprocess.Popen(cmd,
                         shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
            
        # read and print output
        for line in iter(p.stdout.readline, b''):
        #for line in p.stdout.readlines():
            ln = line.decode()
            # line is returned with a \n character at the end 
            # ln = ln[0:len(ln)-2]
            #print(ln, end="", flush=True)
            
        retval = p.wait()
    except:
        pass
    if retval == 0:
        ret = True
        
    return ret
    

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
        qr.printlog("Listing available examples ...", loglevel=1)
    
        import quantarhei.wizard.examples as exmpl
         
        filenames = exmpl._available_examples
        
        if args.glob:
            pattern = args.glob
            matching = _match_filenames(filenames, pattern, add_stars=True)
        else:
            matching = filenames
            
        for ex in matching:
            qr.printlog("    "+ex, loglevel=1)

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
        qr.printlog("Fetching example(s) ...", loglevel=1)
    
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
                    
                    qr.printlog("    "+filename, loglevel=1)

        else:
            qr.printlog("No matching examples found", loglevel=1)

    else:
        parser_fetch.print_help()


def do_command_config(args):
    """Configures Quantarhei
    
    """
    
    qr.printlog("Setting configuration", loglevel=1)


def do_command_report(args):
    """Reports on Quantarhei and the system
    
    """
    
    qr.printlog("Probing system configuration", loglevel=1)

    
def do_command_file(args):
    """Report on the content of a file 
    
    """
    
    qr.printlog("Checking file info", loglevel=1)    
    
    #
    # File name
    # 

    if args.fname:
        fname = args.fname[0]   
        qr.printlog("File name: "+fname, loglevel=1)

    try:
        check = qr.check_parcel(fname)
        qr.printlog("Object type: "+check["class_name"], loglevel=1)
        if check["class_name"] == "builtins.list":
            prcl = qr.load_parcel(fname)
            qr.printlog("List length: "+str(len(prcl)), loglevel=1)
            types = []
            for el in prcl:
                tp = type(prcl[0])
                if tp not in types:
                    types.append(tp)
            qr.printlog("Element types: "+str(types), loglevel=1)
        qr.printlog("Saved with Quantarhei version: "+check["qrversion"],
                    loglevel=1)
        qr.printlog("Description: "+check["comment"])
    except:
        qr.printlog("The file is not a Quantarhei parcel", loglevel=1)
        #print(traceback.format_exc())

        
def do_command_script(args):
    """Provides script management functionality
    
    
    script --path
    
        Reports path to the quantarhei scripts
        
    script --set-path PATH
    
        Sets Quantarhei script path
        
    script --install SCRIPT_FILE
    
        Installs script to the Quantarhei script path
        
    script --list
    
        Lists available scrips
        
    script --history SCRIPT_NAME
    
        Shows the script installation history 
        
    script --status SCRIPT_NAME
    
        Shows script installation status 
        
    script --hash FILE
    
        Reports hash corresponding to the file
        
    
    
        
    
        
    
    """
    global parser_script
    

    pass
    
    
    
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
    parser.add_argument("-y", "--verbosity", type=str, default="5", 
                        help="defines verbosity between 0 and 10")
    parser.add_argument("-f", "--filename", metavar="FILENAME",
                        default="qrhei.log", 
                        help="defines verbosity for logging into file")
    parser.add_argument("-l", "--logtofile", action='store_true', 
                        help="copy logging into a file")
   

 
    
    #
    # Subparser for command `run`
    #
   
    parser_run = subparsers.add_parser("run", help="Script runner")
    
    parser_run.add_argument("script", metavar='script', type=str, 
                          help='script file to be processed', nargs=1)
    parser_run.add_argument("-s", "--silent", action='store_true', 
                          help="logging level set to zero")
    parser_run.add_argument("-q", "--quiet", action='store_true', 
                          help="no output from qrhei script itself")
    parser_run.add_argument("-p", "--parallel", action='store_true', 
                          help="executes the code in parallel")
    parser_run.add_argument("-n", "--nprocesses", type=int, default=0,
                          help="number of processes to start")
    parser_run.add_argument("-f", "--hostfile", metavar="HOSTFILE", 
                            default="", help="list of available"
                            +" host for parallel calculation")
    parser_run.add_argument("-b", "--benchmark", type=int, default=0, 
                          help="run one of the predefined benchmark"
                          +" calculations")
    parser_run.add_argument("-i", "--inputfile", type=str, 
                          default="input.yaml", help="input file for the "
                          +"script")

    
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
    # Subparser for command `file`
    #
    parser_file = subparsers.add_parser("file", help=
                                "Shows information about files"
                                +" created by Quantarhei")
    
    parser_file.add_argument("fname", metavar='fname', type=str, 
                          help='file to be checked', nargs=1) 
    
    parser_file.set_defaults(func=do_command_file)     
    
    
    #
    # Subparser for command `script`
    #

    parser_script = subparsers.add_parser("script", help="Installs and manages"
                                        +" scripts")

    #parser_script.add_argument("glob", metavar='glob', type=str, 
    #                          help='file name', nargs="?")
    #parser_script.add_argument("-i", "--install", action='store_true',
    #                          help="installs a script file")
    parser_script.add_argument("-i", "--install", metavar="SCRIPT", 
                            default="", help="installs a script file") 
    parser_script.add_argument("-s", "--set-path", metavar="PATH", 
                            default="", help="sets path to Quantarhei scripts")    
    parser_script.add_argument("-p", "--path", action='store_true', 
                          help="reports the path to Quantarhei scripts") 
    parser_script.add_argument("-l", "--list", action='store_true', 
                          help="lists all installed scripts")  
    parser_script.add_argument("-y", "--history", metavar="SCRIPT", 
                            default="", help="reports the history of an"+
                            " installed script")
    parser_script.add_argument("-t", "--status", metavar="SCRIPT", 
                            default="", help="reports the status of an"+
                            " installed script")
    parser_script.add_argument("-a", "--hash", metavar="FILE", 
                            default="", help="creates a hash unique to the"
                            +" content of a file")
    parser_script.set_defaults(func=do_command_script)    
     
    
    
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
    
        
    try:      
        if args.func:
            args.func(args)
    except:
        parser.error("No arguments provided")
