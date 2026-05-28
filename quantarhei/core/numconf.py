class NumConf:
    """Configuration of numerical processing in Quantarhei"""

    def __init__(self) -> None:

        # Shared memory acceleration
        self.cpu_acceleration = True
        self.num_threads = -1

        # GPU acceleration
        self.gpu_acceleration = False

        # Libraries
        self.enable_pytorch = False
