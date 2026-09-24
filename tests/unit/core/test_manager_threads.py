"""Thread-local units and basis contexts of the Manager (issue #301).

Each thread must see its own ``energy_units`` and ``eigenbasis_of`` state,
so that calculations running concurrently in different threads cannot
corrupt each other's units or basis stacks.
"""

import threading
import unittest

import numpy

import quantarhei as qr
from quantarhei import Manager, eigenbasis_of, energy_units, set_current_units
from quantarhei.core.units import conversion_facs_energy

from .test_BasisManaged import BasisManagedObject

# conversion factors are exact products of floats; allow a few ulps
RTOL = 1e-12


def _run_in_threads(targets, barrier=None):
    """Run callables in parallel threads and re-raise the first failure."""
    errors = []

    def wrap(fn):
        def run():
            try:
                fn()
            except BaseException as e:  # pragma: no cover - reported below
                errors.append(e)
                if barrier is not None:
                    # release the other threads instead of waiting for timeout
                    barrier.abort()

        return run

    threads = [threading.Thread(target=wrap(fn)) for fn in targets]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=60)
    if errors:
        raise errors[0]
    for t in threads:
        if t.is_alive():  # pragma: no cover
            raise AssertionError("thread did not finish")


class TestManagerThreadLocalContexts(unittest.TestCase):
    def setUp(self):
        set_current_units()
        self.manager = Manager()

    def tearDown(self):
        set_current_units()

    def test_concurrent_energy_units_and_eigenbasis_contexts(self):
        """Two threads hold different unit and basis contexts at the same time"""
        n_rounds = 20
        barrier = threading.Barrier(2, timeout=30)
        manager = self.manager

        h_a = numpy.array([[0.0, 100.0], [100.0, 1000.0]])
        h_b = numpy.array([[0.0, 0.01], [0.01, 0.2]])

        def worker(units, hval, nested):
            cfac = conversion_facs_energy[units]
            for _ in range(n_rounds):
                with energy_units(units):
                    H = qr.Hamiltonian(data=hval)
                    barrier.wait()
                    # the other thread is now inside its own context
                    self.assertEqual(manager.get_current_units("energy"), units)
                    self.assertEqual(manager._in_eu_count, 1)
                    self.assertTrue(manager._in_energy_units_context)
                    numpy.testing.assert_allclose(
                        H._data, hval * cfac, rtol=RTOL, atol=0.0
                    )

                with eigenbasis_of(H):
                    depth = 1
                    if nested:
                        B = qr.Hamiltonian(data=numpy.array([[0.0, 0.3], [0.3, 0.5]]))
                        ctx = eigenbasis_of(B)
                        ctx.__enter__()
                        depth = 2
                    barrier.wait()
                    self.assertEqual(manager.basis_stack, list(range(depth + 1)))
                    self.assertEqual(len(manager.basis_transformations), depth + 1)
                    self.assertTrue(manager._in_eigenbasis_of_context)
                    if not nested:
                        self.assertIs(manager.current_basis_operator, H)
                        # H is diagonal in its own eigenbasis
                        dat = H.data
                        numpy.testing.assert_allclose(
                            dat - numpy.diag(numpy.diag(dat)),
                            0.0,
                            atol=1e-12 * numpy.abs(dat).max(),
                        )
                    barrier.wait()
                    if nested:
                        ctx.__exit__(None, None, None)
                        self.assertIs(manager.current_basis_operator, H)

                barrier.wait()
                # both threads are back outside of all contexts
                self.assertEqual(manager.basis_stack, [0])
                self.assertEqual(manager.basis_registered, {})
                self.assertIsNone(manager.current_basis_operator)
                self.assertFalse(manager._in_eigenbasis_of_context)
                self.assertEqual(manager.get_current_units("energy"), "1/fs")
                # back-transformation from the eigenbasis is exact up to
                # round-off of the orthogonal eigenvector matrix
                numpy.testing.assert_allclose(
                    H._data, hval * cfac, rtol=0.0, atol=1e-12 * abs(hval * cfac).max()
                )

        _run_in_threads(
            [
                lambda: worker("1/cm", h_a, nested=False),
                lambda: worker("eV", h_b, nested=True),
            ],
            barrier=barrier,
        )

    def test_new_thread_starts_outside_of_contexts(self):
        """A new thread does not inherit contexts active in the main thread"""
        seen = {}
        H = BasisManagedObject(numpy.array([[0.1, 1.0], [1.0, 0.0]]), "H")

        def probe():
            m = Manager()
            seen["units"] = m.get_current_units("energy")
            seen["stack"] = list(m.basis_stack)
            seen["eu"] = m._in_energy_units_context
            seen["eb"] = m._in_eigenbasis_of_context
            seen["op"] = m.current_basis_operator

        with energy_units("1/cm"), eigenbasis_of(H):
            _run_in_threads([probe])
            self.assertEqual(self.manager.get_current_units("energy"), "1/cm")
            self.assertEqual(self.manager.basis_stack, [0, 1])

        self.assertEqual(
            seen,
            {"units": "1/fs", "stack": [0], "eu": False, "eb": False, "op": None},
        )

    def test_global_units_seed_new_threads(self):
        """Units set by the module level set_current_units apply to new threads"""
        seen = {}

        def probe():
            seen["units"] = Manager().get_current_units("energy")

        set_current_units({"energy": "1/cm"})
        _run_in_threads([probe])
        self.assertEqual(seen["units"], "1/cm")

        set_current_units()
        _run_in_threads([probe])
        self.assertEqual(seen["units"], "1/fs")

    def test_units_changed_in_thread_do_not_leak(self):
        """Units set inside a worker thread do not change the main thread"""

        def worker():
            m = Manager()
            m.set_current_units("energy", "eV")
            self.assertEqual(m.get_current_units("energy"), "eV")

        _run_in_threads([worker])
        self.assertEqual(self.manager.get_current_units("energy"), "1/fs")


if __name__ == "__main__":
    unittest.main()
