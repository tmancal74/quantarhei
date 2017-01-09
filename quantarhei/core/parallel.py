# -*- coding: utf-8 -*-
"""


"""
class distributed_configuration:

    def __init__(self):
        self.have_mpi = False
        try:
            from mpi4py import MPI
            self.have_mpi = True
        except:
            self.have_mpi = False

        if self.have_mpi:
            comm = MPI.COMM_WORLD
            self.comm = comm
            self.rank = comm.Get_rank()
            self.size = comm.Get_size()
        else:
            self.comm = None
            self.rank = 0
            self.size = 1
            
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
    
if __name__ == "__main__":
 
    N1 = 0
    N2 = 16

    p = distributed_configuration()

    p.print_info()
    for i in block_distributed_range(p, N1, N2):
        print(i, "on rank", p.rank)    
    
