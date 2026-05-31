from __future__ import annotations

import enum


class EnergyUnit(str, enum.Enum):
    INT = "int"
    PER_FS = "1/fs"
    PER_CM = "1/cm"
    EV = "eV"
    MEV = "meV"
    THZ = "THz"
    J = "J"
    SI = "SI"
    NM = "nm"
    HA = "Ha"
    AU = "a.u."


class FrequencyUnit(str, enum.Enum):
    INT = "int"
    PER_FS = "1/fs"
    PER_CM = "1/cm"
    THZ = "THz"
    HZ = "Hz"
    SI = "SI"
    NM = "nm"
    HA = "Ha"
    AU = "a.u."


class LengthUnit(str, enum.Enum):
    INT = "int"
    ANGSTROM = "A"
    NM = "nm"
    BOHR = "Bohr"
    AU = "a.u."
    METER = "m"
    SI = "SI"


class TimeUnit(str, enum.Enum):
    INT = "int"
    FS = "fs"
    AS = "as"
    PS = "ps"
    NS = "ns"
    MS = "ms"
    US = "Ms"
    S = "s"
    SI = "SI"


class TemperatureUnit(str, enum.Enum):
    INT = "int"
    PER_FS = "1/fs"
    KELVIN = "Kelvin"
    CELSIUS = "Celsius"
    PER_CM = "1/cm"
    EV = "eV"
    MEV = "meV"
    THZ = "Thz"
    SI = "SI"
