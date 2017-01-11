# -*- coding: utf-8 -*-
"""


"""
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
            
        self.mp_enabled = True # Massage Passing
        self.host_acc_enabled = True # OpenMP
        self.device_acc_enabled = True # OpenACC, CUDA, or similar
        
        self.manager = Manager()
        
        
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
            print(txt)
            
                
    def reduce(self, A, B, operation="sum"):
        """ Performs a reduction operation on an array
        
        Performs a reduction on an array according to a specification by
        the `operation` argument using MPI. The result is available
        on the process with rank = 0

        """ 
        if operation == "sum":
        
            if not self.have_mpi:
                B += A
                return
            else:
                from mpi4py import MPI
                self.comm.Reduce(A, B, op=MPI.SUM)
                
        else:
            raise Exception("Unknow reduction operation")
            
            
    def allreduce(self, A, B, operation="sum"):
        """ Performs a reduction operation on an array
        
        Performs a reduction on an array according to a specification by
        the `operation` argument using MPI. Result is available to all
        processes.

        """ 

        if operation == "sum":
                       
            if not self.have_mpi:
                B += A
                return
            else:
                from mpi4py import MPI
                self.comm.Allreduce(A, B, op=MPI.SUM)
                
        else:
            raise Exception("Unknow reduction operation")            
            
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
    host_acceleration_enabled = True
    
    if host_acceleration_enabled:
    
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
        
        
if __name__ == "__main__":
 
    import numpy
    
    N1 = 0
    N2 = 16
    
    A = numpy.zeros((N2,N2), dtype=numpy.float64)
    B = numpy.zeros((N2,N2), dtype=numpy.float64)
    dc = DistributedConfiguration()
    dc.print_info()
    for i in block_distributed_range(dc, N1, N2):
        print(i, "on rank", dc.rank)    
        for k in range(N2):
            A[k,k] += i
    dc.reduce(A, B, operation="sum")
    if dc.rank == 0:        
        print(B, ((N2-1)/2)*(N2))
