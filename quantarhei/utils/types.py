from __future__ import annotations

import numbers
from functools import partial
from typing import Any

import numpy


def array_property(name: str, shape: tuple[int, ...] | None = None) -> Any:
    """Controls the access to an array property of package classes


    """
    storage_name = '_'+name

    @property  # type: ignore[misc]
    def prop(self: Any) -> numpy.ndarray:
        return getattr(self,storage_name)

    @prop.setter
    def prop(self: Any, value: Any) -> None:
        try:
            vl = check_numpy_array(value)
            if not (shape == None):
                if not (shape == vl.shape):
                    raise TypeError(
                    f'{name} must be of shape {shape}')
            setattr(self,storage_name,vl)
        except TypeError:
            raise TypeError(
            f'{name} must be either a list or numpy.array')

    return prop


def units_managed_property(name: str, dtype: type) -> Any:
    """Scalar property with units managed

    Warning: The type of the property depends on the object; The object
    has to be EnergyUnitsManaged or similar.
    """
    storage_name = '_'+name

    @property  # type: ignore[misc]
    def prop(self: Any) -> Any:
        val = getattr(self,storage_name)
        return self.convert_2_current_u(val) # This is a method defined in
                                             # the class which handles units
    @prop.setter
    def prop(self: Any, value: Any) -> None:
        if isinstance(value,dtype):
            setattr(self,storage_name,self.convert_2_internal_u(value))
        else:
            raise TypeError(
            '{} must be either a real or complex number')

    return prop



def units_managed_array_property(name: str, dtype: type, shape: tuple[int, ...] | None = None) -> Any:
    """Array property with units managed

    Warning: The type of the property depends on the object; The object
    has to be EnergyUnitsManaged or similar.
    """
    storage_name = '_'+name

    @property  # type: ignore[misc]
    def prop(self: Any) -> Any:
        val = getattr(self,storage_name)
        return self.convert_2_current_u(val) # This is a method defined in
                                             # the class which handles units


    @prop.setter
    def prop(self: Any, value: Any) -> None:

        try:
            vl = check_numpy_array(value)
            if not (shape == None):
                if not (shape == vl.shape):
                    raise TypeError(
                    f'{name} must be of shape {shape}')
            setattr(self,storage_name,self.convert_2_internal_u(vl))
        except TypeError:
            raise TypeError(
            f'{name} must be either a list or numpy.array')

    return prop


def basis_managed_array_property(name: str, dtype: type, shape: tuple[int, ...] | None = None) -> Any:

    storage_name = '_'+name

    @property  # type: ignore[misc]
    def prop(self: Any) -> Any:

        # check if basis id corresponds to the one known by Manager
        cb = self.manager.get_current_basis()
        # get object's current basis
        ob = self.get_current_basis()


        if cb == ob:
            pass
        else:
            # change basis
            self.manager.transform_to_current_basis(self)

        return getattr(self,storage_name)


    @prop.setter
    def prop(self: Any, value: Any) -> None:

        # check if basis id corresponds to the one known by Manager
        cb = self.manager.get_current_basis()
        # get object's current basis
        ob = self.get_current_basis()

        if cb == ob:
            pass
        else:
            # change basis
            self.manager.transform_to_current_basis(self)

        try:
            vl = check_numpy_array(value)
            if not (shape == None):
                if not (shape == vl.shape):
                    raise TypeError(
                    f'{name} must be of shape {shape}')
            setattr(self,storage_name,vl)
        except TypeError:
            raise TypeError(
            f'{name} must be either a list or numpy.array')

    return prop




def managed_array_property(name: str, dtype: type, shape: tuple[int, ...] | None = None) -> Any:
    """Property with the units and basis management


    """
    storage_name = '_'+name

    @property  # type: ignore[misc]
    def prop(self: Any) -> Any:

        # check if basis id corresponds to the one known by Manager
        cb = self.manager.get_current_basis()
        # get object's current basis
        ob = self.get_current_basis()


        if cb == ob:
            pass
        else:
            # change basis
            self.manager.transform_to_current_basis(self)

        val = getattr(self,storage_name)
        return self.convert_2_current_u(val)


    @prop.setter
    def prop(self: Any, value: Any) -> None:

        # check if basis id corresponds to the one known by Manager
        cb = self.manager.get_current_basis()
        # get object's current basis
        ob = self.get_current_basis()


        if cb == ob:
            pass
        else:
            # change basis
            self.manager.transform_to_current_basis(self)

        try:
            vl = check_numpy_array(value)
            if not (shape == None):
                if not (shape == vl.shape):
                    raise TypeError(
                    f'{name} must be of shape {shape}')
            setattr(self,storage_name,self.convert_2_internal_u(vl))
        except TypeError:
            raise TypeError(
            f'{name} must be either a list or numpy.array')

    return prop


def check_numpy_array(val: Any) -> numpy.ndarray:
    """Checks if argument is a numpy array.

    If the argument is a list, it converts it into numpy array. Otherwise,
    error occurs.

    """
    if isinstance(val,numpy.ndarray):
        return val
    if isinstance(val,list):
        try:
            vl = numpy.array(val)
        except TypeError:
            raise TypeError('Numerical array is required')
        return numpy.array(vl)
    raise TypeError('List or numpy.ndarray required')



def typed_property(name: str, dtype: type) -> Any:
    storage_name = '_'+name

    @property  # type: ignore[misc]
    def prop(self: Any) -> Any:
        return getattr(self,storage_name)

    @prop.setter
    def prop(self: Any, value: Any) -> None:
        if isinstance(value,dtype):
            setattr(self,storage_name,value)
        else:
            raise TypeError(f'{name} must be of type {dtype}')

    return prop



def list_of_typed_property(name: str, dtype: type) -> Any:
    storage_name = '_'+name

    @property  # type: ignore[misc]
    def prop(self: Any) -> Any:
        return getattr(self,storage_name)

    @prop.setter
    def prop(self: Any, value: Any) -> None:
        if isinstance(value,list):
            for el in value:
                if not isinstance(el,dtype):
                    raise TypeError(f'{name} must contain \
                    values of type {dtype})',dtype)

    return prop



def alt_type_property(name: str, dtype: Any) -> Any:
    storage_name = '_'+name

    @property  # type: ignore[misc]
    def prop(self: Any) -> Any:
        return getattr(self,storage_name)

    @prop.setter
    def prop(self: Any, value: Any) -> None:
        if isinstance(dtype,list):
            isset = False
            for tps in dtype:
                if isinstance(value,tps):
                    setattr(self,storage_name,value)
                    isset = True
                    break
            if not isset:
                raise TypeError(f'{name} must be of type {dtype}')
        elif isinstance(value,dtype):
            setattr(self,storage_name,value)
        else:
            raise TypeError(f'{name} must be of type {dtype}')

    return prop



def derived_type(name: str, dtype_in: Any) -> Any:
    if isinstance(dtype_in,list):
        tp = partial(alt_type_property,dtype=dtype_in)
    else:
        tp = partial(typed_property,dtype=dtype_in)
    return tp(name)



Float   = partial(typed_property,dtype=numbers.Real)

Integer = partial(typed_property,dtype=numbers.Integral)

Bool    = partial(typed_property,dtype=bool)

BasisManagedComplexArray = partial(basis_managed_array_property,
                                   dtype=numbers.Complex)

BasisManagedRealArray = partial(basis_managed_array_property,
                                dtype=numbers.Real)

UnitsManagedComplexArray = partial(units_managed_array_property,
                                   dtype=numbers.Complex)
UnitsManagedRealArray = partial(units_managed_array_property,
                                dtype=numbers.Real)

UnitsManagedComplex = partial(units_managed_property,
                              dtype=numbers.Complex)

UnitsManagedReal = partial(units_managed_property,
                           dtype=numbers.Real)

ManagedComplexArray = partial(managed_array_property,dtype=numbers.Complex)

ManagedRealArray = partial(managed_array_property,dtype=numbers.Real)
