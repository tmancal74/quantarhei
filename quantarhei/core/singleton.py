from __future__ import annotations

from typing import Any


class Singleton(type):
    """Base type of singletons, such as the main Quantarhei class Manager


    Recipe "Creating a singleton in Python" from Stack Overflow.


    Usage:
    ------
    class MyClass(BaseClass, metaclass=Singleton)

    """

    _instances: dict[type, Singleton] = {}

    def __call__(cls, *args: Any, **kwargs: Any) -> Any:
        if cls not in cls._instances:
            cls._instances[cls] = super().__call__(*args,
                                                                **kwargs)
        return cls._instances[cls]



