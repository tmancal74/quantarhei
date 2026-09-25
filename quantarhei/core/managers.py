"""This class handles several important package wide tasks:

1) Usage of units across objects storing data
2) Basis conversion of all registered objects
3) Calls to proper optimized implementations of numerically heavy
sections of the calculations


Manager is a singleton class, only one instance exists at all times
and all managing objects have the instance of the Manager.

Properies
---------

version : string
contains the package version number


allower_utypes : list
contains a list of unit types which can be controlled by the Manager

units : dictionary
dictionary of available units for each units type

units_repre : dictionary
dictionary of abreviations used to represent various units

units_repre_latex : dictionary
dictionary of latex prepresentations of available units



Units Management
----------------
Units management is performed for all classes derived from
quantarhei.managers.UnitsManaged class.


Basis Conversion Management
---------------------------
Units management is performed for all classes derived from
quantarhei.managers.BasisManaged class.

Basis management works like this: when an class is defined, and its
property needs to be basis managed, one should use a predefined type
`basis_managed_array_property`




"""

from __future__ import annotations

import itertools
import os
import threading
import types
import warnings
import weakref
from abc import ABC, abstractmethod
from collections.abc import Callable
from typing import Any

from ..exceptions import BasisError, ConfigurationError, QuantarheiError, UnitsError


class SecurityWarning(UserWarning):
    """Warning about security-sensitive operations."""


#
# This stops future warnings, notably those in h5py library
# FIXME: remove this in "future"
#
warnings.simplefilter(action="ignore", category=FutureWarning)

import json
from importlib.metadata import PackageNotFoundError as _PackageNotFoundError
from importlib.metadata import version as _pkg_version
from importlib.resources import files

import numpy

from .genconf import GenConf
from .logconf import LogConf
from .numconf import NumConf
from .singleton import Singleton
from .units import (
    conversion_facs_energy,
    conversion_facs_frequency,
    conversion_facs_length,
)

_HARDWIRED_UNITS: dict[str, str] = {
    "energy": "1/fs",
    "frequency": "1/fs",
    "dipolemoment": "Debye",
    "temperature": "Kelvin",
    "length": "A",
}


#
# Basis ids are unique process-wide. Id 0 denotes the default basis, which is
# the same in every thread; every eigenbasis_of context entered in any thread
# receives a fresh id from this counter, so that an object transformed into
# a basis of one thread can never be mistaken for being in a basis of another
# thread. Live ids are mapped to the name of the thread owning the context,
# to produce an informative error message.
#
_basis_ids = itertools.count(1)
_basis_ids_lock = threading.Lock()
_live_bases: dict[int, str] = {}


class _ContextState(threading.local):
    """Per-thread state of the units and basis contexts of the Manager.

    Every thread sees its own instance of these attributes, so that
    ``energy_units`` and ``eigenbasis_of`` contexts entered in one thread do
    not affect calculations running concurrently in another thread. The
    state of a thread is created lazily on its first access; ``units_seed``
    supplies the units a new thread starts with.

    The basis ids on ``basis_stack`` are unique within the process (see
    ``_basis_ids``), so a basis managed object carrying the id of a context
    of another thread is detected and rejected with a :class:`BasisError`
    instead of being silently used in the wrong basis.
    """

    def __init__(self, units_seed: Callable[[], dict[str, str]]) -> None:
        self.current_units: dict[str, str] = units_seed()
        self.saved_units: dict[str, str] = {}
        self.in_energy_units_context = False
        self.in_eu_count = 0
        self.in_eigenbasis_of_context = False
        self.in_eb_count = 0
        self.basis_stack: list[int] = [0]
        self.basis_transformations: list[Any] = [1]
        self.basis_registered: dict[int, weakref.WeakValueDictionary[int, Any]] = {}
        self.current_basis_operator: Any = None


class Manager(metaclass=Singleton):
    """Main package Manager.

    Handles units management, basis conversion, and selection of optimized
    implementations for the entire Quantarhei package. Only one instance
    exists at any time (Singleton pattern).

    The state of the units and basis contexts (``current_units``,
    ``basis_stack``, ``basis_registered``, ``current_basis_operator`` and
    the context flags) is thread-local: each thread has its own stacks, so
    that contexts can be used safely from concurrently running threads.
    Configuration (implementations, ``num_conf`` etc.) is shared by all
    threads.

    A newly started thread begins in the default basis, with the units last
    set by :func:`set_current_units` or :meth:`Manager.set_current_units`
    (internal units unless changed). It does *not* inherit units of an
    ``energy_units`` (or ``frequency_units``, ``length_units``) context
    active in the thread that started it, nor its ``eigenbasis_of``
    contexts; enter the context inside the new thread instead. (Before
    thread-local contexts were introduced, a thread started inside
    ``energy_units("eV")`` saw ``"eV"``.)

    Threads may share basis managed objects (operators, Hamiltonians, ...)
    only while all of them work in the default basis. Using an object that
    another thread has transformed into the basis of its ``eigenbasis_of``
    context raises :class:`~quantarhei.exceptions.BasisError`; basis ids are
    unique within the process, so such an object is never mistaken for
    being in the current thread's basis. The detection is not a lock:
    reading an object while another thread transforms it in place is a data
    race. Give each thread its own copy of the objects it uses inside basis
    contexts (e.g. ``copy.deepcopy``, made while the object is in the
    default basis).

    Pickling or copying the Manager yields the Manager of the current
    process, see :meth:`__reduce__`.

    Attributes
    ----------
    version : str
        The installed package version number.
    allowed_utypes : list of str
        Unit types that can be managed (``'energy'``, ``'frequency'``, etc.).
    units : dict
        Available unit strings for each unit type.
    units_repre : dict
        Short string abbreviations for each unit string.
    units_repre_latex : dict
        LaTeX representations for each unit string.
    """

    try:
        version = _pkg_version("quantarhei")
    except _PackageNotFoundError:
        version = "unknown"

    # hard wired unit options
    allowed_utypes = [
        "energy",
        "frequency",
        "dipolemoment",
        "temperature",
        "time",
        "length",
    ]

    units = {
        "energy": [
            "1/fs",
            "int",
            "1/cm",
            "eV",
            "meV",
            "THz",
            "J",
            "SI",
            "nm",
            "Ha",
            "a.u.",
        ],
        "frequency": ["1/fs", "int", "1/cm", "THz", "Hz", "SI", "nm", "Ha", "a.u."],
        "dipolemoment": ["Debye", "a.u"],
        "temperature": [
            "1/fs",
            "int",
            "Kelvin",
            "Celsius",
            "1/cm",
            "eV",
            "meV",
            "Thz",
            "SI",
        ],
        "time": ["fs", "int", "as", "ps", "ns", "Ms", "ms", "s", "SI"],
        "length": ["int", "A", "nm", "Bohr", "a.u.", "m", "SI"],
    }

    units_repre = {
        "Kelvin": "K",
        "Celsius": "C",
        "Debye": "D",
        "1/cm": "1/cm",
        "THz": "THz",
        "eV": "eV",
        "1/fs": "1/fs",
        "int": "1/fs",
        "meV": "meV",
        "nm": "nm",
        "Ha": "Ha",
        "a.u.": "a.u.",
    }

    units_repre_latex = {
        "Kelvin": "K",
        "Celsius": "C",
        "Debye": "D",
        "1/cm": "cm$^-1$",
        "THz": "THz",
        "eV": "eV",
        "1/fs": "fs$^{-1}$",
        "int": "1/fs",
        "meV": "meV",
        "nm": "nm",
        "Ha": "Ha",
        "a.u.": "a.u.",
    }

    def __reduce__(self) -> tuple[Any, tuple[()]]:
        """Pickle the Manager as a reference to the process-wide singleton

        The Manager holds thread-local state and locks, which cannot be
        pickled, and only one Manager may exist in a process. Unpickling
        (and ``copy.copy``/``copy.deepcopy``) therefore returns the Manager
        of the current process, with its own configuration and the calling
        thread's units and basis contexts; no state of the pickled Manager
        is transferred.
        """
        return (Manager, ())

    def __init__(self) -> None:

        # units with which newly started threads begin
        self._default_units_lock = threading.Lock()
        self._default_units: dict[str, str] = dict(_HARDWIRED_UNITS)

        # thread-local state of units and basis contexts
        self._ctx = _ContextState(self._new_thread_units)

        # main configuration file
        cfile = "~/.quantarhei/quantarhei.json"

        # test the presence of configuration directory
        conf_path = os.path.dirname(cfile)
        self.conf_path = os.path.expanduser(conf_path)
        self.cfile = os.path.expanduser(cfile)

        exists = os.path.exists(self.conf_path)
        isdir = os.path.isdir(self.conf_path)
        if not exists:
            # create directory
            os.mkdir(self.conf_path)

            # write default configuration
            self.main_conf: dict[str, Any] = {
                "units": "units.json",
                "implementations": "implementations.json",
            }

            # save it
            with open(self.cfile, "w") as f:
                json.dump(self.main_conf, f)

        elif exists and (not isdir):
            raise ConfigurationError("Cannot create configuration directory.")

        else:
            # load the main configuration file
            with open(self.cfile) as f:
                self.main_conf = json.load(f)

        #
        # Enforcement of contexts by functions (process-wide switch); the
        # context flags themselves are thread-local, see _ContextState
        #
        self._enforce_contexts = True

        #
        #  Setting physical units
        #

        # internal units are hardwired
        self.internal_units: dict[str, str] = {
            "energy": "1/fs",
            "frequency": "1/fs",
            "dipolemoment": "Debye",
            "temperature": "Kelvin",
            "length": "A",
        }

        # current units are read from conf file
        if not exists:
            # set hard wired defaults and save them
            self.current_units = {
                "energy": "1/fs",
                "frequency": "1/fs",
                "dipolemoment": "Debye",
                "temperature": "Kelvin",
                "length": "A",
            }

            # save them
            self.save_units()

        else:
            self.load_units()

        self.current_units = {
            "energy": "1/fs",
            "frequency": "1/fs",
            "dipolemoment": "Debye",
            "temperature": "Kelvin",
            "length": "A",
        }

        #
        #  Setting implementations
        #

        self.implementation_points: dict[str, str] = {
            "secular-standard-Redfield-rates": "redfield.ssRedfieldRateMatrix"
        }

        #
        #  All available implementations
        #
        self.all_implementations: dict[str, dict[str, str]] = {
            "redfieldrates.ssRedfieldRateMatrix": {
                "0": "quantarhei.implementations.python",
                "1": "quantarhei.implementations.cython",
            }
        }

        self.all_implementations["redfieldtensor.ssRedfieldTensor"] = {
            "0": "quantarhei.implementations.python",
            "1": "quantarhei.implementations.cython",
        }

        self.default_implementations: dict[str, str] = {
            "redfieldrates.ssRedfieldRateMatrix": "0",
            "redfieldtensor.ssRedfieldRateTensor": "0",
        }

        self.optimal_implementations: dict[str, str] = {
            "redfieldrates.ssRedfieldRateMatrix": "1"
        }

        self.current_implementations: dict[str, str] = {
            "redfieldrates.ssRedfieldRateMatrix": "0",
            "redfieldtensor.ssRedfieldRateTensor": "0",
        }

        if not exists:
            # and save them
            self.save_implementations()

        #        else:
        #            self.load_implementations()

        self.change_implementation_at_runtime = True

        self.warn_about_basis_change = False
        self.warn_about_basis_changing_objects = False

        self.save_dict: dict[str, Any] = {}

        #
        # Configuration controlable from qrhei (conf file and qrhei script)
        #
        self.num_conf = NumConf()

        self.log_conf = LogConf()

        self.use_pytorch = False
        self.use_gpu = False

        self.gen_conf = GenConf()

        #
        # Read central configuration from ./quantarhei directory
        #

        #
        # Read local user config file (this will only be done on request)
        #
        # self._read_uconf()

    #
    # Thread-local state of units and basis contexts
    #

    def _new_thread_units(self) -> dict[str, str]:
        """Returns the units a newly started thread begins with"""
        with self._default_units_lock:
            return dict(self._default_units)

    def _set_new_thread_units(self, utype: str, units: str) -> None:
        """Sets the units of type ``utype`` for threads started later"""
        with self._default_units_lock:
            self._default_units[utype] = units

    @property
    def current_units(self) -> dict[str, str]:
        """Units currently used in this thread, per unit type"""
        return self._ctx.current_units

    @current_units.setter
    def current_units(self, value: dict[str, str]) -> None:
        self._ctx.current_units = value

    @property
    def _saved_units(self) -> dict[str, str]:
        return self._ctx.saved_units

    @_saved_units.setter
    def _saved_units(self, value: dict[str, str]) -> None:
        self._ctx.saved_units = value

    @property
    def _in_energy_units_context(self) -> bool:
        return self._ctx.in_energy_units_context

    @_in_energy_units_context.setter
    def _in_energy_units_context(self, value: bool) -> None:
        self._ctx.in_energy_units_context = value

    @property
    def _in_eu_count(self) -> int:
        return self._ctx.in_eu_count

    @_in_eu_count.setter
    def _in_eu_count(self, value: int) -> None:
        self._ctx.in_eu_count = value

    @property
    def _in_eigenbasis_of_context(self) -> bool:
        return self._ctx.in_eigenbasis_of_context

    @_in_eigenbasis_of_context.setter
    def _in_eigenbasis_of_context(self, value: bool) -> None:
        self._ctx.in_eigenbasis_of_context = value

    @property
    def _in_eb_count(self) -> int:
        return self._ctx.in_eb_count

    @_in_eb_count.setter
    def _in_eb_count(self, value: int) -> None:
        self._ctx.in_eb_count = value

    @property
    def basis_stack(self) -> list[int]:
        """Stack of basis ids of this thread; ``0`` is the default basis"""
        return self._ctx.basis_stack

    @basis_stack.setter
    def basis_stack(self, value: list[int]) -> None:
        self._ctx.basis_stack = value

    @property
    def basis_transformations(self) -> list[Any]:
        """Transformation matrices leading to the bases on ``basis_stack``"""
        return self._ctx.basis_transformations

    @basis_transformations.setter
    def basis_transformations(self, value: list[Any]) -> None:
        self._ctx.basis_transformations = value

    @property
    def basis_registered(self) -> dict[int, weakref.WeakValueDictionary[int, Any]]:
        """Weakly held operators to transform back on exit, per basis id"""
        return self._ctx.basis_registered

    @basis_registered.setter
    def basis_registered(
        self, value: dict[int, weakref.WeakValueDictionary[int, Any]]
    ) -> None:
        self._ctx.basis_registered = value

    @property
    def current_basis_operator(self) -> Any:
        """Operator defining the innermost basis context of this thread"""
        return self._ctx.current_basis_operator

    @current_basis_operator.setter
    def current_basis_operator(self, value: Any) -> None:
        self._ctx.current_basis_operator = value

    def load_conf(self) -> None:
        """Loads configuration file

        This is to be called in scripts and notebooks

        """
        self._read_uconf()

    def _read_uconf(self) -> None:
        """Reads user defined local config file

        From Stackoverflow recipe:
            https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path


        """
        fname = self.gen_conf.conf_file_name
        fdir = self.gen_conf.conf_file_path
        fpath = os.path.join(fdir, fname)

        from pathlib import Path

        cfile = Path(fpath)

        if cfile.exists() & cfile.is_file():
            self._load_uconf(fpath)

        else:
            if cfile.exists():
                raise QuantarheiError(
                    "Configuration file "
                    + fpath
                    + " seems to exist"
                    + " but it is not a file"
                )
            else:
                print("Warning: Configuration file " + fpath + " does not exit")
                print("Warning: Placing a default configuration are using it")

                resource_path = "/".join(
                    ("core", "conf", "qrhei.py")
                )  # pragma: no cover
                content = (
                    files("quantarhei").joinpath(resource_path).read_bytes()
                )  # pragma: no cover

                with open(fpath, "w") as f:
                    f.write(content.decode("utf-8"))

                self._load_uconf(fpath)

        # printlog("Configuration file: ", fpath, "loaded", loglevel=9)

    def _load_uconf(self, fpath: str) -> None:
        """ """
        import platform
        import stat

        warnings.warn(
            f"Loading configuration from {fpath} — this file is executed as "
            "Python code. Ensure it comes from a trusted source.",
            SecurityWarning,
            stacklevel=2,
        )

        if platform.system() != "Windows":
            try:
                file_mode = os.stat(fpath).st_mode
                if file_mode & (stat.S_IWGRP | stat.S_IWOTH):
                    warnings.warn(
                        f"Configuration file {fpath} is group- or "
                        "world-writable. Consider running: "
                        f"chmod 600 {fpath}",
                        SecurityWarning,
                        stacklevel=2,
                    )
            except OSError:
                pass

        try:
            import importlib.util

            spec = importlib.util.spec_from_file_location("qrconf", fpath)
            foo = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(foo)
            foo.configure(self)
        except Exception:
            raise QuantarheiError()

    def save_settings(self) -> None:

        # main configuration file
        with open(self.cfile, "w") as f:
            json.dump(self.main_conf, f)

        # units setting
        self.save_units()
        # implementations setting
        self.save_implementations()

    def save_implementations(self) -> None:
        # set the implementations to standard
        implementations = {
            "imp_points": self.implementation_points,
            "all_available": self.all_implementations,
            "default": self.default_implementations,
            "optimal": self.optimal_implementations,
            "current": self.current_implementations,
        }
        imp_file = self.main_conf["implementations"]
        imp_file = os.path.join(self.conf_path, imp_file)
        with open(imp_file, "w") as f:
            json.dump(implementations, f)

    def load_implementations(self) -> None:
        imp_file = self.main_conf["implementations"]
        imp_file = os.path.join(self.conf_path, imp_file)
        with open(imp_file) as f:
            implementations = json.load(f)
            self.implementation_points = implementations["imp_points"]
            self.all_implementations = implementations["all_available"]
            self.default_implementations = implementations["default"]
            self.optimal_implementations = implementations["optimal"]
            self.current_implementations = implementations["current"]

    def save_units(self) -> None:
        units_file = self.main_conf["units"]
        units_file = os.path.join(self.conf_path, units_file)
        with open(units_file, "w") as f:
            json.dump(self.current_units, f)

    def load_units(self) -> None:
        units_file = self.main_conf["units"]
        units_file = os.path.join(self.conf_path, units_file)
        with open(units_file) as f:
            self.current_units = json.load(f)

    def get_real_type(self) -> type:
        """Returns default numpy float type"""
        import numpy

        return numpy.float64

    def get_complex_type(self) -> type:
        """Returns default numpy complex type"""
        import numpy

        return numpy.complex128

    def store_current_basis_operator(self, op: Any) -> None:
        self.current_basis_operator = op

    def remove_current_basis_operator(self) -> None:
        self.current_basis_operator = None

    def unit_repr(self, utype: str = "energy", mode: str = "current") -> str:
        """Returns a string representing the currently used units"""
        if utype in self.allowed_utypes:
            if mode == "current":
                return self.units_repre[self.current_units[utype]]
            if mode == "internal":
                return self.units_repre[self.internal_units[utype]]
            raise UnitsError("Unknown representation mode")

        else:
            raise UnitsError("Unknown unit type")

    def unit_repr_latex(self, utype: str = "energy", mode: str = "current") -> str:
        """Returns a string representing the currently used units"""
        if utype in self.allowed_utypes:
            if mode == "current":
                return self.units_repre_latex[self.current_units[utype]]
            if mode == "internal":
                return self.units_repre_latex[self.internal_units[utype]]
            raise UnitsError("Unknown representation mode")

        else:
            raise UnitsError("Unknown unit type")

    def set_current_units(self, utype: str, units: str) -> None:
        """Set the current units for a given unit type.

        The units are set for the calling thread and, like the module level
        :func:`set_current_units`, also become the starting units of threads
        which access the Manager for the first time afterwards. Threads
        which are already running keep their own units. Use the
        ``energy_units`` and ``length_units`` context managers to change
        units for the calling thread only.

        Parameters
        ----------
        utype : str
            Unit type (e.g. ``'energy'``, ``'length'``).
        units : str
            Unit string to set as current (must be a recognized value
            for the given ``utype``).

        Raises
        ------
        UnitsError
            If ``utype`` is not in ``allowed_utypes`` or ``units`` is not
            recognized for ``utype``.
        """
        self._set_thread_units(utype, units)
        self._set_new_thread_units(utype, units)

    def unset_current_units(self, utype: str) -> None:
        """Restores units saved by the last :meth:`set_current_units`

        The restored units are also used as the starting units of threads
        started afterwards.
        """
        self._unset_thread_units(utype)
        self._set_new_thread_units(utype, self.current_units[utype])

    def _set_thread_units(self, utype: str, units: str, save: bool = True) -> None:
        """Sets the current units of the calling thread only

        With ``save=True`` the previous units are saved and can be restored
        with :meth:`_unset_thread_units`. Units contexts keep their own
        backup and pass ``save=False``, so that a context nested between a
        set/unset pair does not overwrite the saved units.
        """
        if save:
            self._saved_units[utype] = self.get_current_units(utype)

        if utype in self.allowed_utypes:
            if units in self.units[utype]:
                self.current_units[utype] = units
            else:
                raise UnitsError(f"Unknown units of {utype}")
        else:
            raise UnitsError("Unknown type of units")

    def _unset_thread_units(self, utype: str) -> None:
        """Restores previously saved units of the calling thread"""
        try:
            cunits = self._saved_units[utype]
        except KeyError:
            raise UnitsError("Units to restore not found")

        if utype in self.allowed_utypes:
            if cunits in self.units[utype]:
                self.current_units[utype] = cunits
            else:
                raise UnitsError(f"Unknown units of {utype}")
        else:
            raise UnitsError("Unknown type of units")

    def get_current_units(self, utype: str) -> str:
        """ """
        if utype in self.allowed_utypes:
            return self.current_units[utype]
        raise UnitsError("Unknown type of units")

    #    @deprecated
    def cu_energy(  # type: ignore[return]
        self, val: float | numpy.ndarray, units: str = "1/cm"
    ) -> float | numpy.ndarray | None:
        """Converst to current energy units"""
        if units in self.units["energy"]:
            x = conversion_facs_energy[units]
            i_val = x * val

            cu = self.current_units["energy"]
            if cu != "1/fs":
                y = conversion_facs_energy[cu]
                return i_val / y

            return i_val

    #    @deprecated
    def iu_energy(  # type: ignore[return]
        self, val: float | numpy.ndarray, units: str = "1/cm"
    ) -> float | numpy.ndarray | None:
        """Converst to internal energy units"""
        if units in self.units["energy"]:
            x = conversion_facs_energy[units]
            i_val = x * val
            return i_val

    @staticmethod
    def _convert_nm(val: float | numpy.ndarray, cfact: float) -> float | numpy.ndarray:
        tiny = numpy.finfo(float).tiny
        try:
            nonzero = numpy.abs(val) > tiny  # type: ignore[operator]
            ret = numpy.zeros(val.shape, dtype=val.dtype)  # type: ignore[union-attr]
            ret[nonzero] = 1.0 / val[nonzero]  # type: ignore[index]
            return ret / cfact
        except (AttributeError, TypeError):
            return (0.0 if abs(val) <= tiny else 1.0 / val) / cfact  # type: ignore[arg-type]

    def convert_energy_2_internal_u(
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        """Convert energy from currently used units to internal units

        Parameters
        ----------
        val : number, array, list, tuple of numbers
            values to convert

        """
        units = self.current_units["energy"]
        cfact = conversion_facs_energy[units]

        if units == "nm":
            return self._convert_nm(val, cfact)
        return val * cfact

    def convert_energy_2_current_u(
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        """Converts energy from internal units to currently used units

        Parameters
        ----------
        val : number, array, list, tuple of numbers
            values to convert

        """
        units = self.current_units["energy"]
        cfact = conversion_facs_energy[units]

        if units == "nm":
            return self._convert_nm(val, cfact)
        return val / cfact

    def convert_frequency_2_internal_u(
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        """Converts frequency from currently used units to internal units

        Parameters
        ----------
        val : number, array, list, tuple of numbers
            values to convert

        """
        return val * conversion_facs_frequency[self.current_units["frequency"]]

    def convert_frequency_2_current_u(
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        """Converts frequency from internal units to currently used units

        Parameters
        ----------
        val : number, array, list, tuple of numbers
            values to convert

        """
        return val / conversion_facs_frequency[self.current_units["frequency"]]

    def convert_length_2_internal_u(
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        """Converts length from currently used units to internal units

        Parameters
        ----------
        val : number, array, list, tuple of numbers
            values to convert

        """
        return val * conversion_facs_length[self.current_units["length"]]

    def convert_length_2_current_u(
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        """Converts frequency from internal units to currently used units

        Parameters
        ----------
        val : number, array, list, tuple of numbers
            values to convert

        """
        return val / conversion_facs_length[self.current_units["length"]]

    def get_implementation_prefix(self, package: str = "", taskname: str = "") -> str:
        # default_imp_prefix = "quantarhei.implementations.python"

        pname = package + "." + taskname
        whichone = self.current_implementations[pname]
        imp_prefix = self.all_implementations[pname][str(whichone)]

        return imp_prefix

    def get_implementation_points(self) -> dict[str, str]:
        return self.implementation_points

    def get_all_implementations(self) -> dict[str, dict[str, str]]:
        return self.all_implementations

    def get_all_implementations_of(self, imp: str) -> dict[str, str]:
        imp_id = self.implementation_points[imp]
        return self.all_implementations[imp_id]

    def get_current_implementation(self, imp: str) -> str:
        imp_id = self.implementation_points[imp]
        whichone = self.current_implementations[imp_id]
        return self.all_implementations[imp_id][str(whichone)]

    def set_current_implementation(self, imp: str, choice: str) -> None:
        imp_id = self.implementation_points[imp]
        self.current_implementations[imp_id] = choice

    def get_current_basis(self) -> int:
        """Returns the current basis id"""
        return self._ctx.basis_stack[-1]

    def set_new_basis(self, SS: Any) -> int:
        """Pushes a new basis reached by transformation ``SS`` on the stack

        Returns the id of the new basis. Ids are unique within the process
        (not the stack depth), so that ids of different threads never
        coincide.
        """
        ctx = self._ctx
        with _basis_ids_lock:
            nb = next(_basis_ids)
            _live_bases[nb] = threading.current_thread().name
        ctx.basis_stack.append(nb)
        ctx.basis_transformations.append(SS)
        ctx.basis_registered[nb] = weakref.WeakValueDictionary()
        return nb

    def _release_basis(self, bb: int) -> None:
        """Marks the basis ``bb`` as no longer used by any context"""
        with _basis_ids_lock:
            _live_bases.pop(bb, None)

    def transform_to_current_basis(self, operator: Any) -> None:
        """Transforms an operator to the currently used basis

        Parameters
        ----------
        operator : operator
            Any basis managed operator


        """
        ob = operator.get_current_basis()
        cb = self.get_current_basis()

        if self.warn_about_basis_changing_objects:
            print(
                "Object ",
                operator.__class__,
                id(operator),
                " is changing basis from ",
                ob,
                " to: ",
                cb,
            )

        if ob != cb:
            SS = numpy.diag(numpy.ones(operator.dim))
            # find out if current basis of the object is in the stack (i.e. it
            # was used sometime in the past)
            basis_stack = self._ctx.basis_stack
            basis_transformations = self._ctx.basis_transformations
            if ob in basis_stack:
                sl = len(basis_stack)
                # scroll back over the bases
                for k in range(1, sl):
                    # take the basis transformation to the earlier used basis
                    ZZ = basis_transformations[sl - k]

                    # included it into the transformation matrix
                    SS = numpy.dot(ZZ, SS)
                    # if the basis is found, break away from the loop
                    if basis_stack[sl - k - 1] == ob:
                        break
            else:
                with _basis_ids_lock:
                    owner = _live_bases.get(ob)
                if owner is not None:
                    raise BasisError(
                        f"Basis of the object is not on stack: object "
                        f"{operator.__class__.__name__} is held in basis {ob} "
                        f"of an eigenbasis_of context of thread {owner!r}, "
                        f"not of the current thread "
                        f"{threading.current_thread().name!r}. Basis managed "
                        "objects must not be shared between threads while "
                        "any of them is inside an eigenbasis_of context; "
                        "give each thread its own copy of the object, made while "
                        "it is in the default basis."
                    )
                raise BasisError(
                    f"Basis of the object is not on stack: object "
                    f"{operator.__class__.__name__} is held in basis {ob}, "
                    "whose context is no longer active."
                )

            operator.transform(SS)
            operator.set_current_basis(cb)
            self.register_with_basis(cb, operator)

    def register_with_basis(self, nb: int, operator: Any) -> None:
        """Registers an operator to be transformed back when basis ``nb`` exits

        Only a weak reference is kept, so operators which become unreachable
        inside a long-lived context are not kept alive by the Manager.
        Registering an operator which is already registered is a no-op, and
        registration with the default basis ``0`` is ignored, because there
        is no transformation to undo.
        """
        if nb == 0:
            return
        self._ctx.basis_registered[nb][id(operator)] = operator


class Managed:
    """Base class for managed objects"""

    manager = Manager()


class UnitsManaged(Managed):
    """Base class for objects with management of units"""

    def convert_energy_2_internal_u(
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        return self.manager.convert_energy_2_internal_u(val)

    def convert_energy_2_current_u(
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        return self.manager.convert_energy_2_current_u(val)

    def convert_length_2_internal_u(
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        return self.manager.convert_length_2_internal_u(val)

    def convert_length_2_current_u(
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        return self.manager.convert_length_2_current_u(val)

    def unit_repr(self, utype: str = "energy") -> str:
        return self.manager.unit_repr(utype)

    def unit_repr_latex(self, utype: str = "energy") -> str:
        return self.manager.unit_repr_latex(utype)


class _TypedUnitsManaged(Managed):
    _utype: str = ""
    _internal_unit: str = ""

    def convert_2_internal_u(self, val: float | numpy.ndarray) -> float | numpy.ndarray:
        converter = getattr(self.manager, f"convert_{self._utype}_2_internal_u")
        return converter(val)

    def convert_2_current_u(self, val: float | numpy.ndarray) -> float | numpy.ndarray:
        converter = getattr(self.manager, f"convert_{self._utype}_2_current_u")
        return converter(val)

    def unit_repr(self) -> str:
        return self.manager.unit_repr(self._utype)

    def unit_repr_latex(self) -> str:
        return self.manager.unit_repr_latex(self._utype)


class EnergyUnitsManaged(_TypedUnitsManaged):
    _utype = "energy"
    _internal_unit = "1/fs"
    utype = "energy"
    units = "1/fs"


class LengthUnitsManaged(_TypedUnitsManaged):
    _utype = "length"
    _internal_unit = "A"
    utype = "length"
    units = "A"


class BasisManaged(Managed):
    """Base class for objects with managed basis

    The object stores the id of the basis its data are currently expressed
    in. Ids are unique within the process; ``0`` is the default basis shared
    by all threads. An object must not be used by one thread while another
    thread holds it in the basis of its ``eigenbasis_of`` context: this
    raises :class:`~quantarhei.exceptions.BasisError` (see :class:`Manager`).
    """

    _current_basis = Manager().get_current_basis()

    def get_current_basis(self) -> int:
        return self._current_basis

    def set_current_basis(self, bb: int) -> None:
        self._current_basis = bb


class units_context_manager(ABC):
    """General context manager to manage physical units of values"""

    def __init__(self, utype: str = "energy") -> None:
        self.manager = Manager()
        if utype in self.manager.allowed_utypes:
            self.utype = utype
        else:
            raise UnitsError("Unknown units type")

    @abstractmethod
    def __enter__(self) -> None: ...

    @abstractmethod
    def __exit__(
        self,
        ext_ty: type[BaseException] | None,
        exc_val: BaseException | None,
        tb: types.TracebackType | None,
    ) -> None: ...


class energy_units(units_context_manager):
    """Context manager for units of energy.

    Sets the active energy units for the duration of the ``with`` block and
    restores the previous units on exit.

    Parameters
    ----------
    units : str
        Energy unit string (e.g. ``'1/cm'``, ``'eV'``, ``'int'``).

    Examples
    --------
    >>> import quantarhei as qr
    >>> with qr.energy_units("1/cm"):
    ...     H = qr.Hamiltonian(data=[[0.0, 100.0], [100.0, 12000.0]])
    """

    def __init__(self, units: str) -> None:
        super().__init__(utype="energy")

        if units in self.manager.units["energy"]:
            self.units = units
        else:
            raise UnitsError("Unknown energy units")

    def __enter__(self) -> None:
        # save current energy units
        self.units_backup = self.manager.get_current_units("energy")
        self.manager._set_thread_units(self.utype, self.units, save=False)
        self.manager._in_energy_units_context = True
        self.manager._in_eu_count += 1

    def __exit__(
        self,
        ext_ty: type[BaseException] | None,
        exc_val: BaseException | None,
        tb: types.TracebackType | None,
    ) -> None:
        if exc_val is not None and self.units != self.units_backup:
            warnings.warn(
                f"An exception exited a 'with energy_units(\"{self.units}\")' block. "
                f"Unit state has been restored to '{self.units_backup}', but any "
                f"intermediate results computed inside the block may be in unexpected units.",
                stacklevel=2,
            )
        try:
            self.manager._set_thread_units("energy", self.units_backup, save=False)
        finally:
            self.manager._in_eu_count -= 1
            if self.manager._in_eu_count == 0:
                self.manager._in_energy_units_context = False


class frequency_units(energy_units):
    """Context manager for units of frequency.

    Behaves identically to :class:`energy_units` since frequency and energy
    share the same internal representation in Quantarhei.

    Parameters
    ----------
    units : str
        Frequency unit string (e.g. ``'1/cm'``, ``'THz'``, ``'int'``).
    """

    pass


class length_units(units_context_manager):
    """Context manager for length units.

    Sets the active length units for the duration of the ``with`` block and
    restores the previous units on exit.

    Parameters
    ----------
    units : str
        Length unit string (e.g. ``'A'``, ``'nm'``, ``'Bohr'``).
    """

    def __init__(self, units: str) -> None:
        super().__init__(utype="length")

        if units in self.manager.units["length"]:
            self.units = units
        else:
            raise UnitsError("Unknown length units")

    def __enter__(self) -> None:
        # save current energy units
        self.units_backup = self.manager.get_current_units("length")
        self.manager._set_thread_units(self.utype, self.units, save=False)

    def __exit__(
        self,
        ext_ty: type[BaseException] | None,
        exc_val: BaseException | None,
        tb: types.TracebackType | None,
    ) -> None:
        if exc_val is not None and self.units != self.units_backup:
            warnings.warn(
                f"An exception exited a 'with length_units(\"{self.units}\")' block. "
                f"Unit state has been restored to '{self.units_backup}', but any "
                f"intermediate results computed inside the block may be in unexpected units.",
                stacklevel=2,
            )
        self.manager._set_thread_units("length", self.units_backup, save=False)


class basis_context_manager(ABC):
    """General context manager to manage basis"""

    def __init__(self) -> None:
        self.manager = Manager()

    @abstractmethod
    def __enter__(self) -> None: ...

    @abstractmethod
    def __exit__(
        self,
        ext_ty: type[BaseException] | None,
        exc_val: BaseException | None,
        tb: types.TracebackType | None,
    ) -> None: ...


class eigenbasis_of(basis_context_manager):
    """Context manager for working in the eigenbasis of an operator.

    Diagonalizes the given operator on entry, making subsequent calculations
    in the eigenbasis, and transforms all registered objects back on exit.

    Parameters
    ----------
    operator : SelfAdjointOperator
        The operator whose eigenbasis is used inside the context block.
    """

    def __init__(self, operator: Any) -> None:
        super().__init__()
        self.op = operator
        self._previous_basis_operator: Any = None

    def __enter__(self) -> None:

        manager = self.manager

        if manager.warn_about_basis_change:
            print("\nQr >>> Entering basis context manager ...")

        cb = manager.get_current_basis()
        ob = self.op.get_current_basis()

        if cb != ob:
            manager.transform_to_current_basis(self.op)

        # SS = self.op.diagonalize()
        SS = self.op.get_diagonalization_matrix()
        manager.set_new_basis(SS)

        # self.manager.register_with_basis(nb,self.op)
        # self.op.set_current_basis(nb)

        # the operator defining the enclosing context is restored on exit
        self._previous_basis_operator = manager.current_basis_operator
        manager.store_current_basis_operator(self.op)
        manager._in_eigenbasis_of_context = True

        if manager.warn_about_basis_change:
            print("\nQr >>>  ... setting context done")

    def __exit__(
        self,
        ext_ty: type[BaseException] | None,
        exc_val: BaseException | None,
        tb: types.TracebackType | None,
    ) -> None:

        manager = self.manager

        if manager.warn_about_basis_change:
            print("\nQr >>> Returning from basis context manager. Cleaning ...")

        try:
            # This is the basis we are leaving
            bb = manager.basis_stack.pop()
            manager._release_basis(bb)
            # this is the transformation we got here with
            SS = manager.basis_transformations.pop()
            # This is the new basis
            nb = manager.basis_stack[-1]

            # objects registered with the basis we are leaving; the entry is
            # released right away, so that it cannot outlive the context even
            # if the transformation below fails
            registered = manager.basis_registered.pop(bb, None)
            operators = [] if registered is None else list(registered.values())

            # inverse of the transformation matrix
            S1 = numpy.linalg.inv(SS)

            # transform all registered objects
            for op in operators:
                op.transform(S1, inv=SS)
                op.set_current_basis(nb)

                # operators which appeared in this context and were not
                # registered in the one above are now registered with it
                # (registration is idempotent and ignored for basis 0)
                manager.register_with_basis(nb, op)

        finally:
            manager.store_current_basis_operator(self._previous_basis_operator)
            self._previous_basis_operator = None
            if len(manager.basis_stack) == 1:
                manager._in_eigenbasis_of_context = False

        if self.manager.warn_about_basis_change:
            print("\nQr >>> ... cleaning done")


def set_current_units(units: dict[str, str] | None = None) -> None:
    """Set units globally without a context manager.

    The units are set for the calling thread and are also used as the
    starting units of threads which access the Manager for the first time
    afterwards. Threads which are already running keep their own units.
    Units of ``energy_units`` contexts are not passed on to new threads.

    Parameters
    ----------
    units : dict of str, optional
        Mapping from unit type (e.g. ``'energy'``) to unit string
        (e.g. ``'1/cm'``). If ``None``, all unit types are reset to
        their internal defaults.

    Raises
    ------
    Exception
        If any key in ``units`` is not a recognized unit type.
    """
    manager = Manager()
    if units is not None:
        # set units using a supplied dictionary
        for utype in units:
            if utype in manager.allowed_utypes:
                un = units[utype]
                # handle the identity of "frequency" and "energy"
                if utype == "frequency":
                    utype = "energy"
                    un = units["frequency"]

                manager.set_current_units(utype, un)
            else:
                raise UnitsError(f"Unknown units type {utype}")

    else:
        # reset units to the default
        for utype in manager.internal_units:
            if utype in manager.allowed_utypes:
                manager.set_current_units(utype, manager.internal_units[utype])
            else:
                raise UnitsError(f"Unknown units type {utype}")


def units_state() -> dict[str, str]:
    """Returns a snapshot of the current global unit settings.

    Useful for inspecting Manager state in notebooks or debugging
    unexpected unit conversions.

    Returns
    -------
    dict
        Mapping of unit type to current unit string, e.g.
        ``{'energy': '1/cm', 'frequency': '1/fs', 'length': 'A',
        'temperature': 'Kelvin', 'dipolemoment': 'Debye'}``.

    Examples
    --------
    >>> import quantarhei as qr
    >>> state = qr.units_state()
    >>> state['energy']
    '1/fs'

    """
    manager = Manager()
    return dict(manager.current_units)
