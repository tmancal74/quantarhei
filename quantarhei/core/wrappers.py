from __future__ import annotations

import warnings
from collections.abc import Callable
from functools import wraps
from typing import Any

from ..exceptions import QuantarheiError

# import os
from .managers import Manager


def deprecated(
    func: Callable[..., Any] | None = None,
    *,
    alternative: str | None = None,
) -> Callable[..., Any]:
    def decorator(wrapped_func: Callable[..., Any]) -> Callable[..., Any]:
        message = f"{wrapped_func.__qualname__} is deprecated"
        if alternative is not None:
            message = f"{message}; use {alternative} instead"

        @wraps(wrapped_func)
        def wrapper(*arg: Any, **kwargs: Any) -> Any:
            warnings.warn(
                message,
                DeprecationWarning,
                stacklevel=2,
            )
            return wrapped_func(*arg, **kwargs)

        return wrapper

    if func is None:
        return decorator

    return decorator(func)


def prevent_basis_context(func: Callable[..., Any]) -> Callable[..., Any]:
    @wraps(func)
    def wrapper(*arg: Any, **kwargs: Any) -> Any:
        m = Manager()
        if m._in_eigenbasis_of_context and m._enforce_contexts:
            raise QuantarheiError(
                "This function MUST NOT be called"
                " from within an 'eigenbasis_of' context."
            )
        return func(*arg, **kwargs)

    return wrapper


def enforce_basis_context(func: Callable[..., Any]) -> Callable[..., Any]:
    @wraps(func)
    def wrapper(*arg: Any, **kwargs: Any) -> Any:
        m = Manager()
        if (not m._in_eigenbasis_of_context) and m._enforce_contexts:
            raise QuantarheiError(
                "This function MUST be called from within an 'eigenbasis_of' context."
            )
        return func(*arg, **kwargs)

    return wrapper


def prevent_energy_units_context(func: Callable[..., Any]) -> Callable[..., Any]:
    @wraps(func)
    def wrapper(*arg: Any, **kwargs: Any) -> Any:
        m = Manager()
        if m._in_energy_units_context and m._enforce_contexts:
            raise QuantarheiError(
                "This function MUST NOT be called"
                " from within an 'energy_units' context."
            )
        return func(*arg, **kwargs)

    return wrapper


def enforce_energy_units_context(func: Callable[..., Any]) -> Callable[..., Any]:
    @wraps(func)
    def wrapper(*arg: Any, **kwargs: Any) -> Any:
        m = Manager()
        if not m._in_energy_units_context and m._enforce_contexts:
            raise QuantarheiError(
                "This function MUST be called from within an 'energy_units' context."
            )
        return func(*arg, **kwargs)

    return wrapper
