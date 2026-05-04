import numpy
import pytest

import quantarhei as qr
from quantarhei.core.managers import (
    BasisContextError,
    BasisError,
    assert_in_eigenbasis_of,
    assert_not_in_eigenbasis_context,
)


def test_basis_context_error_is_basis_error():
    assert issubclass(BasisContextError, BasisError)


def test_assert_not_in_eigenbasis_context_passes_outside():
    # Should not raise when outside any eigenbasis context
    assert_not_in_eigenbasis_context()


def test_assert_not_in_eigenbasis_context_raises_inside():
    H = qr.Hamiltonian(data=numpy.array([[0.0, 0.1], [0.1, 1.0]]))
    with pytest.raises(BasisContextError):
        with qr.eigenbasis_of(H):
            assert_not_in_eigenbasis_context()


def test_assert_in_eigenbasis_of_raises_outside_context():
    H = qr.Hamiltonian(data=numpy.array([[0.0, 0.1], [0.1, 1.0]]))
    with pytest.raises(BasisError):
        assert_in_eigenbasis_of(H)


def test_assert_in_eigenbasis_of_passes_inside_correct_context():
    H = qr.Hamiltonian(data=numpy.array([[0.0, 0.1], [0.1, 1.0]]))
    with qr.eigenbasis_of(H):
        assert_in_eigenbasis_of(H)  # should not raise


def test_assert_in_eigenbasis_of_raises_for_wrong_operator():
    H1 = qr.Hamiltonian(data=numpy.array([[0.0, 0.1], [0.1, 1.0]]))
    H2 = qr.Hamiltonian(data=numpy.array([[0.0, 0.2], [0.2, 2.0]]))
    with qr.eigenbasis_of(H1):
        with pytest.raises(BasisError):
            assert_in_eigenbasis_of(H2)


def test_data_access_in_correct_basis_does_not_raise():
    H = qr.Hamiltonian(data=numpy.array([[0.0, 0.1], [0.1, 1.0]]))
    # Access in site basis — no context, should work
    _ = H.data


def test_data_access_auto_transforms_in_different_basis():
    H1 = qr.Hamiltonian(data=numpy.array([[0.0, 0.1], [0.1, 1.0]]))
    H2 = qr.Hamiltonian(data=numpy.array([[0.0, 0.2], [0.2, 2.0]]))
    h2_orig = H2.data.copy()
    with qr.eigenbasis_of(H1):
        # H2 auto-transforms to current basis on .data access
        h2_in_eigbasis = H2.data
        # Should be different from original (transformed)
        assert not numpy.allclose(h2_in_eigbasis, h2_orig)
    # Outside, H2 is back to original basis
    assert numpy.allclose(H2.data, h2_orig)


def test_data_set_auto_transforms_in_different_basis():
    H1 = qr.Hamiltonian(data=numpy.array([[0.0, 0.1], [0.1, 1.0]]))
    H2 = qr.Hamiltonian(data=numpy.array([[0.0, 0.2], [0.2, 2.0]]))
    with qr.eigenbasis_of(H1):
        # Setting .data auto-transforms the operator to current basis first
        new_data = numpy.eye(2)
        H2.data = new_data
        # Verify it was set
        assert numpy.allclose(H2.data, new_data)


def test_protect_basis_no_longer_exists():
    H = qr.Hamiltonian(data=numpy.array([[0.0, 0.1], [0.1, 1.0]]))
    assert not hasattr(H, "protect_basis"), (
        "protect_basis should be removed from BasisManaged"
    )
    assert not hasattr(H, "unprotect_basis"), (
        "unprotect_basis should be removed from BasisManaged"
    )
