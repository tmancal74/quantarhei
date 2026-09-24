"""Regression tests for the Manager basis registry (issue #267).

Operators created or transformed inside an ``eigenbasis_of`` context are
registered with the Manager so that they can be transformed back when the
context exits. The registry must not keep otherwise unreachable operators
alive, and it must be released when the context exits.
"""

import gc
import unittest
import weakref

import numpy

from quantarhei import Manager, eigenbasis_of

from .test_BasisManaged import BasisManagedObject


def _live_registered(manager):
    """Number of operators currently held by the registry, per basis id."""
    return {
        bid: len(list(ops.values()) if hasattr(ops, "values") else ops)
        for bid, ops in manager.basis_registered.items()
    }


class TestBasisRegistry(unittest.TestCase):
    def setUp(self):
        self.manager = Manager()
        self.H = BasisManagedObject(numpy.array([[0.1, 1.0], [1.0, 0.0]]), "H")

    def test_registry_empty_outside_contexts(self):
        """Repeated contexts leave no registry entries behind"""
        for _ in range(50):
            with eigenbasis_of(self.H):
                BasisManagedObject(numpy.array([[0.0, 0.3], [0.3, 2.0]]), "A")
        self.assertEqual(self.manager.basis_registered, {})
        self.assertEqual(self.manager.basis_stack, [0])
        self.assertEqual(len(self.manager.basis_transformations), 1)
        self.assertIsNone(self.manager.current_basis_operator)

    def test_discarded_operators_not_accumulated_inside_context(self):
        """Operators dropped inside a long-lived context are not kept alive"""
        n_ops = 200
        with eigenbasis_of(self.H):
            cb = self.manager.get_current_basis()
            before = _live_registered(self.manager).get(cb, 0)
            refs = []
            for _ in range(n_ops):
                op = BasisManagedObject(numpy.array([[0.0, 0.3], [0.3, 2.0]]), "tmp")
                refs.append(weakref.ref(op))
                del op
            gc.collect()
            after = _live_registered(self.manager).get(cb, 0)
            alive = sum(1 for r in refs if r() is not None)

        self.assertEqual(alive, 0)
        self.assertEqual(after, before)

    def test_live_operators_still_transformed_back(self):
        """Weak registration does not drop operators that are still in use"""
        a0 = numpy.array([[0.0, 0.3], [0.3, 2.0]])
        with eigenbasis_of(self.H):
            A = BasisManagedObject(a0.copy(), "A")
            gc.collect()
            with eigenbasis_of(A):
                B = BasisManagedObject(a0.copy(), "B")
                gc.collect()
            gc.collect()
        # A was created in basis 1 with data a0 given in that basis, so outside
        # it equals S a0 S^-1 where S diagonalizes H.
        _, S = numpy.linalg.eigh(self.H.data)
        expected = numpy.dot(S, numpy.dot(a0, numpy.linalg.inv(S)))
        self.assertTrue(numpy.allclose(A.data, expected, rtol=1e-12, atol=1e-12))
        self.assertEqual(A.get_current_basis(), 0)
        self.assertEqual(B.get_current_basis(), 0)

    def test_operator_registered_once_per_basis(self):
        """Nested contexts do not create duplicate registrations"""
        with eigenbasis_of(self.H):
            A = BasisManagedObject(numpy.array([[0.0, 0.3], [0.3, 2.0]]), "A")
            for _ in range(10):
                with eigenbasis_of(A):
                    _ = A.data
            self.assertEqual(_live_registered(self.manager), {1: 1})
        self.assertEqual(self.manager.basis_registered, {})

    def test_registry_released_when_exit_transformation_fails(self):
        """A failing back-transformation does not leave a stale entry"""

        class Broken(BasisManagedObject):
            def transform(self, SS, inv=None):
                if self.get_current_basis() != 0:
                    raise RuntimeError("transform failed")
                super().transform(SS, inv=inv)

        with self.assertRaises(RuntimeError):
            with eigenbasis_of(self.H):
                Broken(numpy.array([[0.0, 0.3], [0.3, 2.0]]), "X")
                keep = Broken(numpy.array([[0.0, 0.3], [0.3, 2.0]]), "Y")

        self.assertEqual(self.manager.basis_registered, {})
        self.assertEqual(self.manager.basis_stack, [0])
        self.assertFalse(self.manager._in_eigenbasis_of_context)
        self.assertIsNone(self.manager.current_basis_operator)

    def test_current_basis_operator_restored_by_nested_exit(self):
        """Leaving an inner context restores the outer basis operator"""
        B = BasisManagedObject(numpy.array([[0.0, 0.3], [0.3, 2.0]]), "B")
        self.assertIsNone(self.manager.current_basis_operator)
        with eigenbasis_of(self.H):
            self.assertIs(self.manager.current_basis_operator, self.H)
            with eigenbasis_of(B):
                self.assertIs(self.manager.current_basis_operator, B)
            self.assertIs(self.manager.current_basis_operator, self.H)
        self.assertIsNone(self.manager.current_basis_operator)

    def test_unentered_context_does_not_hold_operator(self):
        """Constructing eigenbasis_of without entering it has no side effects"""
        ctx = eigenbasis_of(self.H)
        self.assertIsNone(self.manager.current_basis_operator)
        del ctx


if __name__ == "__main__":
    unittest.main()
