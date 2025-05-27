# -*- coding: utf-8 -*-


from quantarhei import Molecule
from quantarhei import energy_units

en = [0.0, 12500] # 12500 1/cm corresponds to an absorption
                  # at the wavelength of 800 nm (near infra red light)
with energy_units("1/cm"):
    m = Molecule(en, name="Molecule with units")


H = m.get_Hamiltonian()
print("\nInternal units:")
print("---------------")
print(H)

with energy_units("1/cm"):
    print("\nInverse centimeters:")
    print("--------------------")
    print(H)

with energy_units("eV"):
    print("\nElectronvolts:")
    print("--------------")
    print("\nManaged data: \n", H.data)
    print("\nRaw data:\n",H._data)    