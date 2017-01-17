# -*- coding: utf-8 -*-
"""


"""

import numpy

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
            
def block_distributed_range(config, start, stop):
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
        
        