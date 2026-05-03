import unittest

import numpy

import quantarhei


class TestTypeConstants(unittest.TestCase):
    def test_REAL_is_float64(self):
        self.assertIs(quantarhei.REAL, numpy.float64)

    def test_COMPLEX_is_complex128(self):
        self.assertIs(quantarhei.COMPLEX, numpy.complex128)

    def test_REAL_is_not_resolved_via_Manager(self):
        self.assertIsNot(type(quantarhei.REAL), type(lambda: None))
        self.assertIs(quantarhei.REAL, numpy.float64)

    def test_constants_usable_as_dtype(self):
        a = numpy.zeros((2, 2), dtype=quantarhei.REAL)
        b = numpy.zeros((2, 2), dtype=quantarhei.COMPLEX)
        self.assertEqual(a.dtype, numpy.float64)
        self.assertEqual(b.dtype, numpy.complex128)
