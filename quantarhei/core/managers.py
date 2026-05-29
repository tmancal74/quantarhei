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

import os
import types
import warnings
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


class Manager(metaclass=Singleton):
    """Main package Manager.

    Handles units management, basis conversion, and selection of optimized
    implementations for the entire Quantarhei package. Only one instance
    exists at any time (Singleton pattern).

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

    def __init__(self) -> None:

        self.current_units: dict[str, str] = {}

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
            self.main_conf: dict[str, Any] = {  # type: ignore[explicit-any]
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

        self.current_basis_operator: Any = None  # type: ignore[explicit-any]

        #
        # Flags for all contexts which are enforced or prevented by functions
        #
        self._enforce_contexts = True
        self._in_eigenbasis_of_context = False
        self._in_eb_count = 0
        self._in_energy_units_context = False
        self._in_eu_count = 0

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

        self.basis_stack: list[int] = []
        self.basis_stack.append(0)
        self.basis_transformations: list[Any] = []  # type: ignore[explicit-any]
        self.basis_transformations.append(1)
        self.basis_registered: dict[int, list[Any]] = {}  # type: ignore[explicit-any]

        self.warn_about_basis_change = False
        self.warn_about_basis_changing_objects = False

        self._saved_units: dict[str, str] = {}

        self.save_dict: dict[str, Any] = {}  # type: ignore[explicit-any]

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

    def store_current_basis_operator(self, op: Any) -> None:  # type: ignore[explicit-any]
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

        Parameters
        ----------
        utype : str
            Unit type (e.g. ``'energy'``, ``'length'``).
        units : str
            Unit string to set as current (must be a recognized value
            for the given ``utype``).

        Raises
        ------
        Exception
            If ``utype`` is not in ``allowed_utypes`` or ``units`` is not
            recognized for ``utype``.
        """
        self._saved_units[utype] = self.get_current_units(utype)

        if utype in self.allowed_utypes:
            if units in self.units[utype]:
                self.current_units[utype] = units
            else:
                raise UnitsError(f"Unknown units of {utype}")
        else:
            raise UnitsError("Unknown type of units")

    def unset_current_units(self, utype: str) -> None:
        """Restores previously saved units of a given type"""
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
    def cu_energy(  # type: ignore[explicit-any, return]
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
    def iu_energy(  # type: ignore[explicit-any, return]
        self, val: float | numpy.ndarray, units: str = "1/cm"
    ) -> float | numpy.ndarray | None:
        """Converst to internal energy units"""
        if units in self.units["energy"]:
            x = conversion_facs_energy[units]
            i_val = x * val
            return i_val

    def convert_energy_2_internal_u(  # type: ignore[explicit-any]
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        """Convert energy from currently used units to internal units

        Parameters
        ----------
        val : number, array, list, tuple of numbers
            values to convert

        """
        units = self.current_units["energy"]
        cfact = conversion_facs_energy[self.current_units["energy"]]

        # special handling for nano meters
        if units == "nm":
            # zero is interpreted as zero energy; use tiny threshold to guard
            # against subnormals that would overflow 1/val to inf
            tiny = numpy.finfo(float).tiny
            try:
                nonzero = numpy.abs(val) > tiny  # type: ignore[operator]
                ret = numpy.zeros(val.shape, dtype=val.dtype)  # type: ignore[union-attr]
                ret[nonzero] = 1.0 / val[nonzero]  # type: ignore[index]
                return ret / cfact
            except (AttributeError, TypeError):
                return (0.0 if abs(val) <= tiny else 1.0 / val) / cfact  # type: ignore[arg-type]
        else:
            return val * cfact

    def convert_energy_2_current_u(  # type: ignore[explicit-any]
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

        # special handling for nanometers
        if units == "nm":
            # zero is interpreted as zero energy; use tiny threshold to guard
            # against subnormals that would overflow 1/val to inf
            tiny = numpy.finfo(float).tiny
            try:
                nonzero = numpy.abs(val) > tiny  # type: ignore[operator]
                ret = numpy.zeros(val.shape, dtype=val.dtype)  # type: ignore[union-attr]
                ret[nonzero] = 1.0 / val[nonzero]  # type: ignore[index]
                return ret / cfact
            except (AttributeError, TypeError):
                return (0.0 if abs(val) <= tiny else 1.0 / val) / cfact  # type: ignore[arg-type]
        else:
            return val / cfact

    def convert_frequency_2_internal_u(  # type: ignore[explicit-any]
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        """Converts frequency from currently used units to internal units

        Parameters
        ----------
        val : number, array, list, tuple of numbers
            values to convert

        """
        return val * conversion_facs_frequency[self.current_units["frequency"]]

    def convert_frequency_2_current_u(  # type: ignore[explicit-any]
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        """Converts frequency from internal units to currently used units

        Parameters
        ----------
        val : number, array, list, tuple of numbers
            values to convert

        """
        return val / conversion_facs_frequency[self.current_units["frequency"]]

    def convert_length_2_internal_u(  # type: ignore[explicit-any]
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        """Converts length from currently used units to internal units

        Parameters
        ----------
        val : number, array, list, tuple of numbers
            values to convert

        """
        return val * conversion_facs_length[self.current_units["length"]]

    def convert_length_2_current_u(  # type: ignore[explicit-any]
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

    def register_implementation(  # type: ignore[explicit-any]
        self, imp_point: str, prefix: str, asint: Any = None
    ) -> None:
        pass

    def commit_implementation(  # type: ignore[explicit-any]
        self, imp_point: str, prefix: str, asint: Any = None
    ) -> None:
        pass

    def get_current_basis(self) -> int:
        """Returns the current basis id"""
        l = len(self.basis_stack)
        return self.basis_stack[l - 1]

    def set_new_basis(self, SS: Any) -> int:  # type: ignore[explicit-any]
        nb = self.get_current_basis() + 1
        self.basis_stack.append(nb)
        self.basis_transformations.append(SS)
        self.basis_registered[nb] = []
        return nb

    def transform_to_current_basis(self, operator: Any) -> None:  # type: ignore[explicit-any]
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
            if ob in self.basis_stack:
                sl = len(self.basis_stack)
                # scroll back over the bases
                for k in range(1, sl):
                    # take the basis transformation to the earlier used basis
                    ZZ = self.basis_transformations[sl - k]

                    # included it into the transformation matrix
                    SS = numpy.dot(ZZ, SS)
                    # if the basis is found, break away from the loop
                    if self.basis_stack[sl - k - 1] == ob:
                        break
            else:
                raise BasisError("Basis of the object is not on stack.")

            operator.transform(SS)
            operator.set_current_basis(cb)
            self.register_with_basis(cb, operator)

    def register_with_basis(self, nb: int, operator: Any) -> None:  # type: ignore[explicit-any]
        self.basis_registered[nb].append(operator)


class Managed:
    """Base class for managed objects"""

    manager = Manager()


class UnitsManaged(Managed):
    """Base class for objects with management of units"""

    def convert_energy_2_internal_u(  # type: ignore[explicit-any]
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        return self.manager.convert_energy_2_internal_u(val)

    def convert_energy_2_current_u(  # type: ignore[explicit-any]
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        return self.manager.convert_energy_2_current_u(val)

    def convert_length_2_internal_u(  # type: ignore[explicit-any]
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        return self.manager.convert_length_2_internal_u(val)

    def convert_length_2_current_u(  # type: ignore[explicit-any]
        self, val: float | numpy.ndarray
    ) -> float | numpy.ndarray:
        return self.manager.convert_length_2_current_u(val)

    def unit_repr(self, utype: str = "energy") -> str:
        return self.manager.unit_repr(utype)

    def unit_repr_latex(self, utype: str = "energy") -> str:
        return self.manager.unit_repr_latex(utype)


class EnergyUnitsManaged(Managed):
    utype = "energy"
    units = "1/fs"

    def convert_2_internal_u(self, val: float | numpy.ndarray) -> float | numpy.ndarray:  # type: ignore[explicit-any]
        return self.manager.convert_energy_2_internal_u(val)

    def convert_2_current_u(self, val: float | numpy.ndarray) -> float | numpy.ndarray:  # type: ignore[explicit-any]
        return self.manager.convert_energy_2_current_u(val)

    def unit_repr(self) -> str:
        return self.manager.unit_repr("energy")

    def unit_repr_latex(self, utype: str = "energy") -> str:
        return self.manager.unit_repr_latex(utype)


class LengthUnitsManaged(Managed):
    """Class providing functions for length units conversion"""

    utype = "length"
    units = "A"

    def convert_2_internal_u(self, val: float | numpy.ndarray) -> float | numpy.ndarray:  # type: ignore[explicit-any]
        return self.manager.convert_length_2_internal_u(val)

    def convert_2_current_u(self, val: float | numpy.ndarray) -> float | numpy.ndarray:  # type: ignore[explicit-any]
        return self.manager.convert_length_2_current_u(val)

    def unit_repr(self) -> str:
        return self.manager.unit_repr(self.utype)

    def unit_repr_latex(self) -> str:
        return self.manager.unit_repr_latex(self.utype)


class BasisManaged(Managed):
    """Base class for objects with managed basis"""

    _current_basis = Manager().get_current_basis()

    def get_current_basis(self) -> int:
        return self._current_basis

    def set_current_basis(self, bb: int) -> None:
        self._current_basis = bb


class units_context_manager:
    """General context manager to manage physical units of values"""

    def __init__(self, utype: str = "energy") -> None:
        self.manager = Manager()
        if utype in self.manager.allowed_utypes:
            self.utype = utype
        else:
            raise UnitsError("Unknown units type")

    def __enter__(self) -> None:
        pass

    def __exit__(
        self,
        ext_ty: type[BaseException] | None,
        exc_val: BaseException | None,
        tb: types.TracebackType | None,
    ) -> None:
        pass


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
        self.manager.set_current_units(self.utype, self.units)
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
            self.manager.set_current_units("energy", self.units_backup)
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
        self.manager.set_current_units(self.utype, self.units)

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
        self.manager.set_current_units("length", self.units_backup)


class basis_context_manager:
    """General context manager to manage basis"""

    def __init__(self) -> None:
        self.manager = Manager()

    def __enter__(self) -> None:
        pass

    def __exit__(
        self,
        ext_ty: type[BaseException] | None,
        exc_val: BaseException | None,
        tb: types.TracebackType | None,
    ) -> None:
        pass


class eigenbasis_of(basis_context_manager):
    """Context manager for working in the eigenbasis of an operator.

    Diagonalizes the given operator on entry, making subsequent calculations
    in the eigenbasis, and transforms all registered objects back on exit.

    Parameters
    ----------
    operator : SelfAdjointOperator
        The operator whose eigenbasis is used inside the context block.
    """

    def __init__(self, operator: Any) -> None:  # type: ignore[explicit-any]
        super().__init__()
        self.op = operator
        self.manager.store_current_basis_operator(self.op)

    def __enter__(self) -> None:

        self.manager._in_eigenbasis_of_context = True

        if self.manager.warn_about_basis_change:
            print("\nQr >>> Entering basis context manager ...")

        cb = self.manager.get_current_basis()
        ob = self.op.get_current_basis()

        if cb != ob:
            self.manager.transform_to_current_basis(self.op)

        # SS = self.op.diagonalize()
        SS = self.op.get_diagonalization_matrix()
        self.manager.set_new_basis(SS)

        # self.manager.register_with_basis(nb,self.op)
        # self.op.set_current_basis(nb)

        if self.manager.warn_about_basis_change:
            print("\nQr >>>  ... setting context done")

    def __exit__(
        self,
        ext_ty: type[BaseException] | None,
        exc_val: BaseException | None,
        tb: types.TracebackType | None,
    ) -> None:

        if self.manager.warn_about_basis_change:
            print("\nQr >>> Returning from basis context manager. Cleaning ...")

        try:
            # This is the basis we are leaving
            bb = self.manager.basis_stack.pop()
            # this is the transformation we got here with
            SS = self.manager.basis_transformations.pop()
            # This is the new basis
            bss = len(self.manager.basis_stack)
            nb = self.manager.basis_stack[bss - 1]

            # inverse of the transformation matrix
            S1 = numpy.linalg.inv(SS)

            # transform all registered objects
            operators = self.manager.basis_registered[bb]

            if nb != 0:
                # operators registered with the context above this one
                ops_above = self.manager.basis_registered[nb]

            for op in operators:
                op.transform(S1, inv=SS)
                op.set_current_basis(nb)

                # operators which appeared in this context and where not
                # register in the one above are now registerd
                if nb != 0:
                    if op not in ops_above:
                        self.manager.register_with_basis(nb, op)

            self.manager.remove_current_basis_operator()

            del self.manager.basis_registered[bb]
        finally:
            if len(self.manager.basis_stack) == 1:
                self.manager._in_eigenbasis_of_context = False

        if self.manager.warn_about_basis_change:
            print("\nQr >>> ... cleaning done")


def set_current_units(units: dict[str, str] | None = None) -> None:
    """Set units globally without a context manager.

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
