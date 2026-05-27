import os

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy

numpy.set_printoptions(precision=8, sign=" ", legacy="1.13")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow-running")
