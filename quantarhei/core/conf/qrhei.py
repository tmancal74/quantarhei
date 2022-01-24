"""
###############################################################################

       Quantarhei computation and logging configuration file
       
###############################################################################

Purpose of this file
--------------------

This file, if placed to the working directory of the computation (default 
search path), is read, and its ``configure`` function is called the
first time ``Manager`` class is instantiated. 


Name and location of Quantarhei configuration file
--------------------------------------------------

Location of the configuration file is controlled by the ``conf_file_path``
attribute of the ``gen_conf`` object of the Manager.

If the file is not found in the prescribed path, default file is placed
to that location instead, and default setting are used (the file can be editted
before the next run). The name of the file is controlled by the
``conf_file_name`` attribute of the ``gen_conf`` object of the Manager
(default is ``qrhei.conf``)


About Quantarhei configuration
------------------------------

Configuration of Quantarhei is hierarchical. Firts, hard defaults are
set in the Manager constructor. Then  configuration from the ``.quantarhei``
directory are read. These override the hard default. Next the ``qrhei.conf``
or its equivalent is read. The settings are overridden again. As the last point
of configuration before running the script, the command line options of the
``qrhei`` driver script are used, again overridding all previous methods of
configuratio.

When the script starts, some configuration parameters cannot be
changed any more (or more precisely, they can be changed, but take no effect).
However, many options can be changed during the execution of the script from
the script iself. For instance, one can change settings for the number of MPI
processes, but as they are set at the start up, this change takes no effect.
On the other hand, GPU usage can be switched on and off during the
computation.

"""

def configure(manager):
    """Configuration of quantarhei computation and logging
    
    """

    ###########################################################################
    # configuration of numerics
    ###########################################################################
    conf = manager.num_conf # DO NOT EDIT THIS LINE
    ###########################################################################

    # usage of mpi
    conf.mpi_acceleration = False 
    
    # here one can prevent competing multithreating if necessary
    conf.cpu_acceleration = True
    conf.num_threads = -1 # -1 means to be determined optimally
        
    # use gpu acceleration (if gpu available)
    # this requires pytorch, but it does
    # not switch on pytorch usage for
    # other than GPU computations
    conf.gpu_acceleration = False 
    
    # restrict yourself only to certain GPUs
    conf.available_gpus = [0, 1]
    
    # enables pytorch as an alternative
    # to numpy even without GPUs
    conf.enable_pytorch = False 
    

    ###########################################################################
    # logging configuration
    ###########################################################################
    conf = manager.log_conf # DO NOT EDIT THIS LINE
    ###########################################################################
    conf.log_on_screen = True
    conf.log_to_file = False
    #conf.log_file_name = "./qrhei.log"
    
    # verbosity is a number from 0 to 10
    # 0 == no information written
    # 10 == all information is written
    conf.verbosity = 5 
    conf.verbose=True

    ###########################################################################
    # general configuration
    ###########################################################################
    conf = manager.gen_conf # DO NOT EDIT THIS LINE
    ###########################################################################
    conf.legacy_relaxation = False
