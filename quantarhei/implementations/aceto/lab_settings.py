from __future__ import annotations

import numpy


class lab_settings:
    FOUR_WAVE_MIXING = 4

    def __init__(self, exptype: int) -> None:
        self.exptype = exptype
        self.orient_aver: numpy.ndarray | None = None

    def set_laser_polarizations(
        self, e1: numpy.ndarray, e2: numpy.ndarray, e3: numpy.ndarray, e4: numpy.ndarray
    ) -> None:
        self.p1 = e1
        self.p2 = e2
        self.p3 = e3
        self.p4 = e4
        # initialiying M4
        M4 = numpy.array(
            [[4.0, -1.0, -1.0], [-1.0, 4.0, -1.0], [-1.0, -1.0, 4.0]],
            dtype=numpy.float64,
        )
        M4 /= 30.0
        F4 = numpy.zeros(3)
        F4[0] = numpy.dot(e4, e3) * numpy.dot(e2, e1)
        F4[1] = numpy.dot(e4, e2) * numpy.dot(e3, e1)
        F4[2] = numpy.dot(e4, e1) * numpy.dot(e3, e2)
        self.orient_aver = numpy.dot(numpy.transpose(M4), F4)

    def oafactor(
        self, d1: numpy.ndarray, d2: numpy.ndarray, d3: numpy.ndarray, d4: numpy.ndarray
    ) -> numpy.ndarray:
        F4 = numpy.zeros(3, dtype=numpy.float64)
        F4[0] = numpy.dot(d4, d3) * numpy.dot(d2, d1)
        F4[1] = numpy.dot(d4, d2) * numpy.dot(d3, d1)
        F4[2] = numpy.dot(d4, d1) * numpy.dot(d3, d2)
        return numpy.dot(self.orient_aver, F4)
