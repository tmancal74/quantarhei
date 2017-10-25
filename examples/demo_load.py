# -*- coding: utf-8 -*-

from quantarhei import load
from quantarhei import read_info

from quantarhei import Molecule


m = Molecule(elenergies=[0.0, 1.0])
m.name = "My molecule"

m.position = [0.0, 1.0, 0.0]


m.save("molecule1.hdf5")


m2 = load("molecule1.hdf5")

info = read_info("molecule1.hdf5")

print(info)