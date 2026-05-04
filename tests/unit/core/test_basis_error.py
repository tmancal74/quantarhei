import numpy
import pytest

import quantarhei as qr
from quantarhei.core.managers import (
    BasisContextError,
    BasisError,
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
