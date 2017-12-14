# -*- coding: utf-8 -*-

class NumConf:
    """Configuration of numerical processing in Quantarhei
    
    """
    
    def __init__(self):
        
        # Message Passing Parallelization
        self.mpi_acceleration = False
        
        # Shared memory parallelization
        self.cpu_acceleration = True
        self.num_threads = -1
        
        # GPU acceleration
        self.gpu_acceleration = False
        
        # Libraries
        self.enable_pytorch = False
        