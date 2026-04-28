from __future__ import annotations

from collections.abc import Callable
from functools import wraps
from importlib import import_module
from typing import Any

from .managers import Manager


def implementation(
    package: str = "",
    taskname: str = "",
    at_runtime: bool = False,
    fallback_local: bool = False,
    always_local: bool = False,
) -> Callable[[Any], Any]:
    """Decorator to select numerical implememtation"""
    m = Manager()

    def decorate_at_runtime(func: Callable[..., Any]) -> Callable[..., Any]:
        """Decoration at run time

        The wrapper decides which function to return at runtime.

        """

        @wraps(func)
        def wrapper(*arg: Any, **kwargs: Any) -> Any:
            fc = get_function(
                func,
                package,
                taskname,
                default_local=fallback_local,
                always_local=always_local,
            )
            return fc(*arg, **kwargs)

        return wrapper

    def decorate_at_loadtime(func: Callable[..., Any]) -> Callable[..., Any]:
        """Decoration at load time

        The wrapper decides which function to return when the Manager module
        is loaded, i.e. at the start of the application.

        """
        fc = get_function(
            func,
            package,
            taskname,
            default_local=fallback_local,
            always_local=always_local,
        )

        @wraps(func)
        def wrapper(*arg: Any, **kwargs: Any) -> Any:
            return fc(*arg, **kwargs)

        return wrapper

    if at_runtime and m.change_implementation_at_runtime:
        return decorate_at_runtime

    return decorate_at_loadtime


#
#  Auxiliary function
#


def load_function(lib: str, fce: str) -> Callable[..., Any]:
    """Load the module and get the desired function"""
    try:
        a = import_module(lib)
    except ImportError:
        print("Cannot load module", lib)

    if hasattr(a, fce):
        fc = getattr(a, fce)
    else:
        raise Exception(f"Cannot reach implementation of {fce} ")

    return fc


def get_function(
    func: Callable[..., Any],
    package: str,
    taskname: str,
    default_local: bool,
    always_local: bool,
) -> Callable[..., Any]:
    """Decide which function to use"""
    if always_local:
        return func

    m = Manager()
    # default implementation package
    default_imp_prefix = "quantarhei.implementations.python"

    # decide which implementation will be used
    imp_prefix = m.get_implementation_prefix(package=package, taskname=taskname)

    # load the package
    try:
        imp_name = imp_prefix + "." + package
        fc = load_function(imp_name, taskname)

    except Exception:
        try:
            # fall back on pure Python implementation
            if default_local:
                fc = func
            else:
                imp_name = default_imp_prefix + "." + package
                fc = load_function(imp_name, taskname)

            # FIXME: issue a warning
            print("WARNING: import failed, falling back on pure Python")
        except Exception:
            # do not provide implementation, call the decorated function itself
            # FIXME: issue a warning (this is an unwanted result)
            print("WARNING: calling decorated function itself")
            fc = func

    return fc
