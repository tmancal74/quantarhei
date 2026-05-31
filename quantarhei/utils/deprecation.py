from __future__ import annotations

import functools
import warnings
from collections.abc import Callable
from typing import Any, TypeVar

F = TypeVar("F", bound=Callable[..., Any])


def deprecated(message: str) -> Callable[[F], F]:
    def decorator(func: F) -> F:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            warnings.warn(
                f"{func.__qualname__} is deprecated and will be removed in a "
                f"future version. {message}",
                DeprecationWarning,
                stacklevel=2,
            )
            return func(*args, **kwargs)

        return wrapper  # type: ignore[return-value]

    return decorator
