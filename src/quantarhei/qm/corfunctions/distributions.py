from __future__ import annotations

import numpy


class BoseEinsteinDistribution:
    def __init__(self, freq_axis: object, temperature: float) -> None:
        kBT = temperature
        self.data = 1.0 / (numpy.exp(freq_axis.data / kBT) - 1.0)
