# -*- coding: utf-8 -*-
"""


"""
import sys
import numpy

from .. import COMPLEX


def call_finish():
    pass


class DistributedConfiguration:

    def __init__(self, schedule=None):

        self.use_steerer = False
        if self.use_steerer:
            nsteerer = 1
        else:
            nsteerer = 0
            
            
        self.have_mpi = False
        try:
            from mpi4py import MPI
            self.have_mpi = True
        except:
            self.have_mpi = False
            #print("WARNING: MPI import failed - switching to serial.")
            #print("WARNING: If you run this with mpiexec/mpirun," 
            #+ " you get independent copies of the serial program")

            
        if self.have_mpi:
            comm = MPI.COMM_WORLD
            self.comm = comm
            self.rank = comm.Get_rank()
            self.stearer = comm.Get_size()
            self.size = self.stearer-nsteerer
        else:
            self.comm = None
            self.rank = 0
            self.size = 1
                    
        self.parallel_level = 0
        
        #if self.size > 1:
        #    self.parallel_level = 1
            
        self.inparallel = False
        self.parallel_region = 0
        
        self.silent = True
        
        
    def stearer_loop(self):
        """ Here we wait for stearer tasks
        
        """
        for ii in range(self.size):
            s = self.comm.recv()
            if s == "FINISH":
                count += 1
                sys.exit()
        else:
            pass
        
        
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
            
        #
        # The process with largest rank enters a stearer loop where it
        # will wait message until it comes
        #
        if self.use_steerer:
            if self.rank == self.size:
                self.steerer_loop()
            else:
                import atexit
                atexit.register(call_finish)
                
        # this counting is independent of whether we have MPI
        self.parallel_region += 1
        

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
            
        self.parallel_region -= 1
        
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

        if self.parallel_region < 1:
            raise Exception("This code has to be run from a declared parallel_region")
 
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
        
        if self.parallel_region < 1:
            raise Exception("This code has to be run from a declared parallel_region")
        
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
        if self.parallel_region < 1:
            raise Exception("This code has to be run from a declared parallel_region")
        
        return self.comm.bcast(value, root=root)

            
def start_parallel_region():
    """Starts a parallel region
    
    This is a clean solution without having to explicitely invoke Manager
    and DistributedConfiguration class
    
    """
    from .managers import Manager
    dc = Manager().get_DistributedConfiguration()
    dc.start_parallel_region()
    if dc.rank != 0:
        Manager().log_conf.verbosity -= 2
        Manager().log_conf.fverbosity -= 2


def close_parallel_region():
    """Closes a parallel region
    
    This is a clean solution without having to explicitely invoke Manager
    and DistributedConfiguration class
    
    """     
    from .managers import Manager
    dc = Manager().get_DistributedConfiguration()
    dc.finish_parallel_region()
    if dc.rank != 0:
        Manager().log_conf.verbosity += 2
        Manager().log_conf.fverbosity += 2

    
def distributed_configuration():
    """Returns the DistributedConfiguration object from the Manager
    
    """
    from .managers import Manager
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
    
    from .managers import Manager
    # we share the work only in parallel_level == 1  
    config = Manager().get_DistributedConfiguration()

    if config.parallel_region < 1:
        raise Exception("This code has to be run from a declared parallel_region")
        
    if config.parallel_level==1:
        
        config.inparallel_entered = True
        
        rng = _calculate_ranges(config, start, stop)
        config.range = rng

        return range(rng[0],rng[1])
        
    else:

        return range(start, stop)


def block_distributed_list(dlist, return_index=False):
    """ Creates sublists for each process 
        
        Returns an iterator over a part of the list
        Altogether the iterators span the original list.

        Parameters
        ----------
        
        dlist : list
            the list to loop over
            
        return_index : bool
            if True, the function returns a list of tuples, where the first
            element is an index counting the interated list elements from 0

    """
    
    from .managers import Manager
    # we share the work only in parallel_level == 1  
    config = Manager().get_DistributedConfiguration()
 
    if config.parallel_region < 1:
        raise Exception("This code has to be run from a declared parallel_region")
 
    if config.parallel_level==1:

        config.inparallel_entered = True
        
        rng = _calculate_ranges_list(config, dlist)
    
        config.range = rng

        if return_index:
            lst = []
            for a in range(rng[0],rng[1]):
                lst.append((a, dlist[a]))
            return lst
        else:
            return dlist[rng[0]:rng[1]]
        
    else:

        rng = [0, len(dlist)]
        config.range = rng
        
        if return_index:
            lst = []
            for a in range(len(dlist)):
                lst.append((a, dlist[a]))
            return lst            
        else: 
            return dlist


def block_distributed_array(array, return_index=False):
    """
    
    """
    
    from .managers import Manager
    # we share the work only in parallel_level == 1  
    config = Manager().get_DistributedConfiguration()

    if config.parallel_region < 1:
        raise Exception("This code has to be run from a declared parallel_region")
    
    if config.parallel_level==1:

        config.inparallel_entered = True
        
        rng = _calculate_ranges_array(config, array)
    
        config.range = rng
        
        if return_index:
            lst = []
            for a in range(array.shape[0]):
                lst.append((a, array[a]))
            return lst             
        else:
            return array[rng[0]:rng[1]]
        
    else:

        rng = [0, array.shape[0]]
        config.range = rng
        
        if return_index:
            lst = []
            for a in range(array.shape[0]):
                lst.append((a, array[a]))
            return lst            
        else: 
            return array


def collect_block_distributed_data(containers, setter_function,
                                  retriever_function, tags=None):
    """Collects distributed data into a container on rank 0 nod
    
    Collects the "data" properties of the objects of the 
    data_class type into a container on the rank = 0 nod.
       
       
    Parameters
    ----------
       
    containers : list
        A list of two container. Firts container is a new empty container
        to which everything will be collected (on nod 0), the other
        is the local container with the data to be set to nod 0
        
    tag_type : {int, str}
        Type of the tag identifying data
        
    setter_function : function
        Function that sets the data to the locally created class
        based on `data` and `tag` recieved from other nods
    
    retriever_function : function
        Function which retrieves the data from the container based on a tag
    
        
    """

    from .managers import Manager
    # we share the work only in parallel_level == 1  
    config = Manager().get_DistributedConfiguration()

    if config.parallel_region < 1:
        raise Exception("This code has to be run from a declared parallel_region")
    
    if config.parallel_level==1:
     
        
        if config.rank == 0:
            
            # collect all data and tags into dictionaries

            rng = config.range
            
            data_shape = (1,1)
            data_type = COMPLEX
                            
            for ii in range(config.size):                     
                rng = config.ranges[ii]
                #print("recieving from:", ii)
                for a in range(rng[0],rng[1]):
                    
                    if tags is not None:
                        tag = tags[a]
                    else:
                        tag = a
                        
                    if ii == 0:
                        # retrieving locally calculated data
                        data = retriever_function(containers[1], tag)
                        # get the type and dimension for later use
                        data_shape = data.shape
                        data_type = data.dtype
                    else:
                        # recieving remotely calculated data
                        data = numpy.zeros(data_shape, dtype=data_type)
                        #data[0] = float(a)
                        #print("recieving", a, "from", ii)
                        config.comm.Recv(data, source=ii, tag=a)
                        
                    #print("setting", a, data)
                    setter_function(containers[0], tag, data)
                
        else:
            
            # send the data to nod 0
            rng = config.range
            #print(config.rank, "sends to nod 0:", rng)

            for a in range(rng[0],rng[1]):
                # sending locally calculated data
                if tags is not None:
                    tag = tags[a]
                else:
                    tag = a
                data = retriever_function(containers[1], tag)
                #print("sending", a, "from", config.rank,"to 0")
                config.comm.Send(data, dest=0, tag=a)
                

    else:
        # if there is no parallelization, we do nothing because all data
        # is already in the contaner
        rng = config.range
        for a in range(rng[0], rng[1]):
            if tags is not None:
                tag = tags[a]
            else:
                tag = a
            data = retriever_function(containers[1], tag)
            setter_function(containers[0], tag, data)
        

def _calculate_ranges(config, start, stop):
    """Calculate which part of a give range should belong to which process


    Parameters
    ----------
    
    config : DistributedConfiguration 
        object holding information about the parallel environment of quantarhei
   
    start : int
        start of the integer range
        
    stop : int
        end of the integer range
    
    """    

    whole_range = stop-start
    per_worker = whole_range // config.size
    remainder = whole_range % config.size  
    
    ranges = [None]*config.size
    for rank in range(config.size):
        N1_local = rank*per_worker
        N2_local = N1_local+per_worker
    
        if rank <= remainder:
            if rank != 0:
                N1_local += rank-1
                N2_local += rank
        else:
            N1_local += remainder
            N2_local += remainder
                
        rng = list()
        rng.append(N1_local)
        rng.append(N2_local)
        ranges[rank] = rng
        
    config.ranges = ranges
    return ranges[config.rank]

    
def _calculate_ranges_list(config, dlist):
    """Calculate which part of a given list should belong to which process

    Parameters
    ----------
    
    config : DistributedConfiguration 
        object holding information about the parallel environment of quantarhei

    dlist : list
        list which will be distributed
        
    """     
    ln = len(dlist)
    start = 0
    stop = ln
    return _calculate_ranges(config, start, stop)


def _calculate_ranges_array(config, array):
    """Calculate which part of a given array should belong to which process

    Parameters
    ----------
    
    config : DistributedConfiguration 
        object holding information about the parallel environment of quantarhei

    array : numpy.array
        an array which will be distributed
        
    """     
    ln = array.shape[0]
    start = 0
    stop = ln
    return _calculate_ranges(config, start, stop)
      

import threading

_sentinel = object()

def _send_to_other_by_MPI(config, i, dest, finish=False):
    
    to_send = numpy.zeros(2, dtype=int)
    to_send[0] = i
    if finish:
        to_send[1] = 1
    else:
        to_send[1] = 0
        
    # send it
    config.comm.Send(to_send, dest=dest)


def _receive_from_MPI(config, source):
    
    # receive it
    data = numpy.zeros(2, dtype=int)
    config.comm.Recv(data, source=source)
    
    if data[1] == 1:
        return None
    else:
        return data[0]



class RangeDistributor(threading.Thread):
    
    def __init__(self, start, stop, queue, config):
        threading.Thread.__init__(self)
        self.threadID = 0
        self.name = "Range Distributor"
        self.counter = 0
        self.strt = start
        self.stop = stop
        
        self.queue = queue
        self.config = config
    
    def run(self):
        
        nproc = self.config.size
        
        l = 0
        for val in range(self.strt, self.stop):
            
            #print("In loop:", val)
            
            if l == 0:
                # send to rank 0 via a queue
                #print("putting", val, "into queue")
                self.queue.put(val)
            else:
                # send to others via MPI
                #print("sending", val, "by MPI")
                _send_to_other_by_MPI(self.config, val, l)
            
            l += 1
            if l >= nproc:
                l = 0
                
        for k in range(nproc):
            
            if l == 0:
                #print("Putting finish signal to queue")
                # send finish signal to rank 0 via a queue
                self.queue.put(_sentinel)
            else:
                #print("Sending finish signal by MPI")
                # send finish signal to others via MPI
                _send_to_other_by_MPI(self.config, 0, l, finish=True)
            
            l += 1
            if l >= nproc:
                l = 0           
        
        
def _get_next_from_thread(queue):

    val = queue.get()
        
    if val is _sentinel:
        
        return None

    return val

    
def _get_next_from_MPI(config, source):
    
    # receive value through MPI
    return _receive_from_MPI(config, source)



def asynchronous_range(start, stop):
    """Range distributing numbers asynchronously among processes
    
    """
    from .managers import Manager
    config = Manager().get_DistributedConfiguration()
    
    if config.parallel_region < 1:
        raise Exception("This code has to be run from a declared parallel_region")    
    
    if config.parallel_level==1:  
        
        from queue import Queue
        q = Queue()
        
        if config.rank == 0:
            
            # start a separate thread which will devide job between processes
            # including the present one
            rd = RangeDistributor(start, stop, q, config)
            rd.start()
            
        not_done = True
        while not_done:
            
            if config.rank == 0:
                
                next_val = _get_next_from_thread(q)
            
            else:
                
                next_val = _get_next_from_MPI(config, 0)
                
            #print("Got: ", next_val, "in", config.rank)
            if next_val is not None:
                #print("Yielding", next_val, "in", config.rank)
                yield next_val
        
            else:
                not_done = False
        
        if config.rank == 0:
            
            # at the end of the work, the Thread rejoins the main process
            rd.join()
            
            
    
    else:

        # serial implementation of the same
        for aa in range(start, stop):
            yield aa
    



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

#    These examples, when tested, will be run in serial mode and the test will
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
        
        
        from .managers import Manager
        
        # FIXME: do we have to be parallel_level aware????
        self.fce = function
        self.leader = leader_rank
        
        manager = Manager()
        self.dc = manager.get_DistributedConfiguration()
        if self.dc.parallel_region < 1:
            raise Exception("This code has to be run from a declared parallel_region")
    
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
            
        
