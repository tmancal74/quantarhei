"""Thread-local units and basis contexts of the Manager (issues #267, #301).

Each thread must see its own ``energy_units`` and ``eigenbasis_of`` state,
so that calculations running concurrently in different threads cannot
corrupt each other's units or basis stacks. Basis ids are unique within the
process, so that a basis managed object shared between threads is never
silently used in the basis of another thread's context: this must raise a
BasisError instead.
"""

import threading
import unittest

import numpy

import quantarhei as qr
from quantarhei import (
    Manager,
    eigenbasis_of,
    energy_units,
    length_units,
    set_current_units,
)
from quantarhei.core.units import conversion_facs_energy
from quantarhei.exceptions import BasisError

from .test_BasisManaged import BasisManagedObject

# conversion factors are exact products of floats; allow a few ulps
RTOL = 1e-12

# basis transformations of 2x2 matrices with O(1) entries by orthogonal
# eigenvector matrices: round-off is a few ulps of the largest entry
ATOL = 1e-12

SHARED_MSG = "must not be shared between threads"


def _eigvecs(a):
    """Eigenvector matrix as computed by the test objects (numpy.linalg.eigh)"""
    return numpy.linalg.eigh(a)[1]


def _in_basis(x, ss):
    """Matrix ``x`` expressed in the basis reached by transformation ``ss``"""
    return numpy.linalg.inv(ss) @ x @ ss


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

        stacks = {}

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
                    stack = list(manager.basis_stack)
                    stacks.setdefault(units, []).append(stack)
                    self.assertEqual(stack[0], 0)
                    self.assertEqual(len(stack), depth + 1)
                    self.assertEqual(len(set(stack)), depth + 1)
                    self.assertEqual(manager.get_current_basis(), stack[-1])
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

        # basis ids are unique process-wide: apart from the default basis 0,
        # no id is ever used by both threads
        ids = {
            units: {b for stack in st for b in stack[1:]}
            for units, st in stacks.items()
        }
        self.assertEqual(ids["1/cm"] & ids["eV"], set())
        self.assertEqual(len(ids["1/cm"]), n_rounds)
        self.assertEqual(len(ids["eV"]), 2 * n_rounds)

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
            self.assertEqual(len(self.manager.basis_stack), 2)

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
        """Units set inside a worker thread do not change running threads"""
        seen = {}

        def worker():
            m = Manager()
            m.set_current_units("energy", "eV")
            self.assertEqual(m.get_current_units("energy"), "eV")

        def probe():
            seen["units"] = Manager().get_current_units("energy")

        _run_in_threads([worker])
        self.assertEqual(self.manager.get_current_units("energy"), "1/fs")
        # like the module level function, the method sets the global units
        _run_in_threads([probe])
        self.assertEqual(seen["units"], "eV")

    def test_manager_method_seeds_new_threads(self):
        """Manager.set_current_units and unset_current_units seed new threads"""
        seen = {}

        def probe():
            m = Manager()
            seen["energy"] = m.get_current_units("energy")
            seen["length"] = m.get_current_units("length")

        self.manager.set_current_units("energy", "1/cm")
        self.manager.set_current_units("length", "nm")
        _run_in_threads([probe])
        self.assertEqual(seen, {"energy": "1/cm", "length": "nm"})

        self.manager.unset_current_units("energy")
        self.manager.unset_current_units("length")
        self.assertEqual(self.manager.get_current_units("energy"), "1/fs")
        _run_in_threads([probe])
        self.assertEqual(seen, {"energy": "1/fs", "length": "A"})

    def test_units_contexts_do_not_seed_new_threads(self):
        """Units of energy_units/length_units contexts stay in their thread

        This is a deliberate change: before the contexts became thread-local,
        a thread started inside ``energy_units("eV")`` saw ``"eV"``.
        """
        seen = {}

        def probe():
            m = Manager()
            seen["energy"] = m.get_current_units("energy")
            seen["length"] = m.get_current_units("length")

        set_current_units({"energy": "1/cm"})
        with energy_units("eV"), length_units("nm"):
            _run_in_threads([probe])
            self.assertEqual(seen, {"energy": "1/cm", "length": "A"})
        _run_in_threads([probe])
        self.assertEqual(seen, {"energy": "1/cm", "length": "A"})

    def test_aggregate_build_does_not_change_global_units(self):
        """Temporary internal units used by Aggregate.build stay thread-local"""
        seen = {}

        def probe():
            seen["energy"] = Manager().get_current_units("energy")

        with energy_units("1/cm"):
            mols = [qr.Molecule([0.0, 12000.0]), qr.Molecule([0.0, 12100.0])]
            agg = qr.Aggregate(molecules=mols)
            agg.set_resonance_coupling(0, 1, 100.0)
            agg.build()
            # the units saved by build are not overwritten by units contexts
            # entered during the build
            self.assertEqual(self.manager.get_current_units("energy"), "1/cm")
        _run_in_threads([probe])
        self.assertEqual(seen["energy"], "1/fs")
        self.assertEqual(self.manager.get_current_units("energy"), "1/fs")


class TestSharedObjectsAcrossThreads(unittest.TestCase):
    """Basis managed objects shared by threads inside basis contexts

    An object transformed into the basis of one thread's ``eigenbasis_of``
    context must never be used as if it were in another thread's basis. The
    guaranteed behaviour is a BasisError in the other thread, while the
    owning thread keeps getting correct data.
    """

    h0 = numpy.array([[0.1, 1.0], [1.0, 0.0]])
    h1 = numpy.array([[0.0, 0.3], [0.3, 2.0]])
    h2 = numpy.array([[1.0, -0.7], [-0.7, 0.2]])
    x0 = numpy.array([[0.5, 0.2], [0.2, -1.0]])

    def setUp(self):
        set_current_units()
        self.manager = Manager()

    def tearDown(self):
        set_current_units()

    def test_shared_context_operator(self):
        """Two threads using eigenbasis_of the same shared H"""
        H = BasisManagedObject(self.h0.copy(), "H")
        a_inside = threading.Event()
        b_done = threading.Event()

        def owner():
            with eigenbasis_of(H):
                # H is now transformed into this thread's eigenbasis
                numpy.testing.assert_allclose(
                    H.data, numpy.diag(numpy.linalg.eigvalsh(self.h0)), atol=ATOL
                )
                a_inside.set()
                self.assertTrue(b_done.wait(30))
                # the failed attempts of the other thread did not touch H
                numpy.testing.assert_allclose(
                    H.data, numpy.diag(numpy.linalg.eigvalsh(self.h0)), atol=ATOL
                )

        def other():
            try:
                self.assertTrue(a_inside.wait(30))
                m = Manager()
                with self.assertRaisesRegex(BasisError, SHARED_MSG):
                    with eigenbasis_of(H):
                        pass  # pragma: no cover
                # the failed __enter__ left no context behind
                self.assertEqual(m.basis_stack, [0])
                self.assertIsNone(m.current_basis_operator)
                self.assertFalse(m._in_eigenbasis_of_context)
                with self.assertRaisesRegex(BasisError, SHARED_MSG):
                    _ = H.data
            finally:
                b_done.set()

        _run_in_threads([owner, other])

        # once the owner left its context, H can be used by any thread again
        numpy.testing.assert_allclose(H.data, self.h0, atol=ATOL)

        def reuse():
            with eigenbasis_of(H):
                numpy.testing.assert_allclose(
                    H.data, numpy.diag(numpy.linalg.eigvalsh(self.h0)), atol=ATOL
                )

        _run_in_threads([reuse])
        numpy.testing.assert_allclose(H.data, self.h0, atol=ATOL)

    def test_shared_operator_in_different_eigenbases(self):
        """A shared X read by threads in eigenbases of different Hamiltonians

        Before basis ids were unique both contexts had id 1 and the second
        thread silently received X in the first thread's basis.
        """
        X = BasisManagedObject(self.x0.copy(), "X")
        a_read = threading.Event()
        b_done = threading.Event()

        def first():
            H1 = BasisManagedObject(self.h1.copy(), "H1")
            with eigenbasis_of(H1):
                numpy.testing.assert_allclose(
                    X.data, _in_basis(self.x0, _eigvecs(self.h1)), atol=ATOL
                )
                a_read.set()
                self.assertTrue(b_done.wait(30))
                numpy.testing.assert_allclose(
                    X.data, _in_basis(self.x0, _eigvecs(self.h1)), atol=ATOL
                )

        def second():
            try:
                H2 = BasisManagedObject(self.h2.copy(), "H2")
                with eigenbasis_of(H2):
                    self.assertTrue(a_read.wait(30))
                    with self.assertRaisesRegex(BasisError, SHARED_MSG):
                        _ = X.data
                    # objects of this thread are unaffected
                    numpy.testing.assert_allclose(
                        H2.data,
                        numpy.diag(numpy.linalg.eigvalsh(self.h2)),
                        atol=ATOL,
                    )
            finally:
                b_done.set()

        _run_in_threads([first, second])
        numpy.testing.assert_allclose(X.data, self.x0, atol=ATOL)

    def test_nested_contexts_in_two_threads(self):
        """Nested contexts in two threads: own objects correct, shared raise"""
        X = BasisManagedObject(self.x0.copy(), "X")
        barrier = threading.Barrier(2, timeout=30)
        a_read = threading.Event()
        b_done = threading.Event()
        stacks = {}

        def worker(name, ha, hb, reads_shared):
            m = Manager()
            HA = BasisManagedObject(ha.copy(), "HA")
            HB = BasisManagedObject(hb.copy(), "HB")
            Y = BasisManagedObject(self.x0.copy(), "Y")
            # transformations the Manager composes for the nested contexts:
            # SB diagonalizes HB as expressed in the eigenbasis of HA
            sa = _eigvecs(ha)
            sb = _eigvecs(_in_basis(hb, sa))
            try:
                with eigenbasis_of(HA):
                    with eigenbasis_of(HB):
                        stacks[name] = list(m.basis_stack)
                        barrier.wait()
                        numpy.testing.assert_allclose(
                            Y.data, _in_basis(self.x0, sa @ sb), atol=ATOL
                        )
                        if reads_shared:
                            numpy.testing.assert_allclose(
                                X.data, _in_basis(self.x0, sa @ sb), atol=ATOL
                            )
                            a_read.set()
                            self.assertTrue(b_done.wait(30))
                        else:
                            self.assertTrue(a_read.wait(30))
                            try:
                                with self.assertRaisesRegex(BasisError, SHARED_MSG):
                                    _ = X.data
                            finally:
                                b_done.set()
                    # back in the eigenbasis of HA only
                    numpy.testing.assert_allclose(
                        Y.data, _in_basis(self.x0, sa), atol=ATOL
                    )
                    if reads_shared:
                        numpy.testing.assert_allclose(
                            X.data, _in_basis(self.x0, sa), atol=ATOL
                        )
                numpy.testing.assert_allclose(Y.data, self.x0, atol=ATOL)
                self.assertEqual(m.basis_stack, [0])
                self.assertEqual(m.basis_registered, {})
            finally:
                a_read.set()
                b_done.set()

        _run_in_threads(
            [
                lambda: worker("a", self.h0, self.h1, True),
                lambda: worker("b", self.h1, self.h2, False),
            ],
            barrier=barrier,
        )
        self.assertEqual(len(stacks["a"]), 3)
        self.assertEqual(len(stacks["b"]), 3)
        self.assertEqual(set(stacks["a"]) & set(stacks["b"]), {0})
        numpy.testing.assert_allclose(X.data, self.x0, atol=ATOL)

    def test_error_after_context_exit(self):
        """An unregistered object left in an exited basis raises clearly"""
        import copy

        H = BasisManagedObject(self.h0.copy(), "H")
        with eigenbasis_of(H):
            _ = H.data  # transforms H into the eigenbasis
            # copies are not registered with the context
            C = copy.deepcopy(H)
        with self.assertRaisesRegex(BasisError, "no longer active"):
            _ = C.data


if __name__ == "__main__":
    unittest.main()
