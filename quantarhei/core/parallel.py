# -*- coding: utf-8 -*-
"""


"""

import numpy

from .managers import Manager

class DistributedConfiguration:

    def __init__(self, schedule=None):

        self.have_mpi = False
        try:
            from mpi4py import MPI
            self.have_mpi = True
        except:
            self.have_mpi = False
            #print("WARNING: MPI import failed - switching to serial.")
            #print("WARNING: If you run this with mpiexec," 
            #+ " you get independent copies of the serial program")

            
        if self.have_mpi:
            comm = MPI.COMM_WORLD
            self.comm = comm
            self.rank = comm.Get_rank()
            self.size = comm.Get_size()
        else:
            self.comm = None
            self.rank = 0
            self.size = 1
                    
        self.parallel_level = 0
        self.inparallel = False
        
        self.silent = True
        
    def start_parallel_region(self):
        """Starts a parallel region
        
        This function raises parallel_level if MPI is present and it sets
        inparallel atribute to True. The distributed ranges provided for 
        parallel for loops are evaluated in parallel only if in parallel_level
        equal 1. This means that tasks are parallelized only from the toppest
        most level. if parallel region is encountered again, no redistribution
        of the tasks is done.
        
        """
        if self.have_mpi:
            if self.size > 1:
                self.parallel_level += 1
            
        if self.parallel_level > 0:
            self.inparallel = True

    def finish_parallel_region(self):
        """Closes a parallel region
        
        Lowers parallel_level by 1
        
        """
        if self.have_mpi:
            if self.size > 1:
                if self.parallel_level == 1:
                    self.comm.Barrier()
                self.parallel_level -= 1
                
        if self.parallel_level < 0:
            raise Exception()
            
        
    def print_info(self):
        if self.rank == 0:
            if self.have_mpi:
                from mpi4py import MPI
                print("Distributed Environemt Info")
                print("------------------------")
                print("Size of the environment : ", self.size)
                print("Processor name          : ", MPI.Get_processor_name())
            else:
                print("No Distributed Environemt")
                
    def print(self, txt, who=0):
        if self.rank == who:
            if not self.silent:
                print(txt)
            
                
    def reduce(self, A, operation="sum"):
        """ Performs a reduction operation on an array
        
        Performs a reduction on an array according to a specification by
        the `operation` argument using MPI. The result is available
        on the process with rank = 0

        """ 
        # only in parallel_level == 1 we share the work
        if self.parallel_level != 1:
            return A

        if operation == "sum":
        
            from mpi4py import MPI
            B = numpy.zeros(A.shape, dtype=A.dtype)
            self.comm.Reduce(A, B, op=MPI.SUM)
            return B    
                
        else:
            raise Exception("Unknown reduction operation")
            
            
    def allreduce(self, A, operation="sum"):
        """ Performs a reduction operation on an array
        
        Performs a reduction on an array according to a specification by
        the `operation` argument using MPI. Result is available to all
        processes.

        """ 

        # only in parallel_level == 1 we share the work
        if self.parallel_level != 1:
            return 
        
        if operation == "sum":
                       
            from mpi4py import MPI
            B = numpy.zeros(A.shape, dtype=A.dtype)
            self.comm.Allreduce(A, B, op=MPI.SUM)
            A[:,:] = B
            
        else:
            raise Exception("Unknown reduction operation")      
            
    def bcast(self, value, root=0):
        #if self.parallel_level != 1:
        #    return value
            
        return self.comm.bcast(value, root=root)

            
def start_parallel_region():
    """Starts a parallel region
    
    This is a clean solution without having to explicitely invoke Manager
    and DistributedConfiguration class
    
    """
    dc = Manager().get_DistributedConfiguration()
    dc.start_parallel_region()


def close_parallel_region():
    """Closes a parallel region
    
    This is a clean solution without having to explicitely invoke Manager
    and DistributedConfiguration class
    
    """    
    dc = Manager().get_DistributedConfiguration()
    dc.finish_parallel_region()
    
    
def distributed_configuration():
    """Returns the DistributedConfiguration object from the Manager
    
    """
    return Manager().get_DistributedConfiguration()

    
def block_distributed_range(start, stop):
    """ Creates an iterator which returns a block of indices
        
        Returns an iterator over a block of indices for each parallel process.
        Altogether the iterators span the range from start to stop.

        Parameters
        ----------
        
        start : int
            beginning of the range
            
        stop: int
            end of the range

    """

    # we share the work only in parallel_level == 1  
    config = Manager().get_DistributedConfiguration()
    if config.parallel_level==1:

        config.inparallel_entered = True
        whole_range = stop-start
        per_worker = whole_range // config.size
        remainder = whole_range % config.size
    
        N1_local = config.rank*per_worker
        N2_local = N1_local+per_worker
        
        if config.rank <= remainder:
            if config.rank != 0:
                N1_local += config.rank-1
                N2_local += config.rank
        else:
            N1_local += remainder
            N2_local += remainder
        
        rng = list()
        rng.append(N1_local)
        rng.append(N2_local)
    
        config.range = rng

        return range(rng[0],rng[1])
        
    else:

        return range(start, stop)

def _parallel_function_wrapper(func, root):
    # FIXME: return a wrapped function with parameter broadcasting
    
    def retfce(params):

        wait_for_work = True
        dc = distributed_configuration()
        wait_for_work = dc.bcast(wait_for_work, root=root)
        params = dc.bcast(params, root=root)
        
        return func(params)
        
    return retfce
    
        
class parallel_function:
    """Context manager handling division of work on parallel tasks 
    
    
    This context manager is motivated by the wish to have a code as
    clean of explicit handling of parallelism as possible. 
    
    
    Parameters
    ----------
    
    function : callable 
        function which will be repeatedly called within the context
    
    leader_rank : int
        rank of the process which handles all the inteligent work
    
        
    Examples
    --------

#    These examples, when tested, will be run serial mode and the test will
#    only confirm that the serial version works
#    
#    
#    Here is some function that uses parallel execution to obtain
#    a single value
#    
#    >>> def pfunction(parameters):
#    ...    val = get_some_value_using_all_processes(parameters)
#    ...    return val
#    
#    
#    And this is a function which takes a function and runs it
#    repetitively examining its return value and deciding when to exit
#    based on the returned value (or some other criterion)
#        
#    >>> def some_inteligence_at_work(func):
#    ...    
#    ...    cont = True
#    ...    val = initial_val
#    ...    while cont:
#    ...        par = select_new_parameters(val)
#    ...        val = func(par)
#    ...        cont = shall_we_continue(val)
#    ...    last_val = val    
#    ...        
#    ...    return last_val
#            
#    
#    Only the leader process will execute the inteligent function and
#    get its return value (this might or might not be a problem)
#           
#    >>> val = 0
#    ...
#    ... with parallel_function(pfunction, leader_rank=0) as f, execute_block:
#    ...     if execute_block:
#    ...         val = some_inteligence_at_work(f)
#
#    
#    If it is a problem that the value val is assigned a final value 
#    on the leader process only, we have to explicitely broadcast it
#    AFTER the context is left
#    
#    >>> dc = Manager().get_DistributedConfiguration()            
#    >>> val = dc.bcast(val, from=0)
        
    """    

    
    def __init__(self, function, leader_rank=0, parallel_level=0):
        
        # FIXME: do we have to be parallel_level aware????
        self.fce = function
        self.leader = leader_rank
        
        manager = Manager()
        self.dc = manager.get_DistributedConfiguration()

    
    def __enter__(self):
        """All except of the leader are put on hold
        
        All processes will wait for the leader (rank=0 by default) 
        to broadcast a list of function arguments.
        
        """
        # helper processes will wait here and should not execute the block
        if self.dc.rank != self.leader:
            
            wait_for_work = True
            # keep waiting for the work in a loop
            while wait_for_work:
                params = None
                # first we communicate if work will come at all
                wait_for_work = self.dc.comm.bcast(wait_for_work,
                                                   root=self.leader)
                if wait_for_work:
                    # now we get parameters for the function to which we bypass
                    params = self.dc.comm.bcast(params, root=self.leader)
                    # call the function with the broadcasted parameters
                    res = self.fce(params)
                    # returned results are irrelevant, everything is handled
                    # by the leader process

            execute_block = False
            ftoret = None
            
        # the leader process goes straight through and will execute the block            
        else:
            
            execute_block = True
            
            ftoret = _parallel_function_wrapper(self.fce, root=self.leader)

        return (ftoret, execute_block)  
    
    def __exit__(self,ext_ty,exc_val,tb):
        """On exit the leader process signals to stop waiting for more work
        
        """
        
        # all processes except of the leader go straight through
        if self.dc.rank == self.leader:
            wait_for_work = False
            wait_for_work = self.dc.bcast(wait_for_work,
                                          root=self.leader)
            
        