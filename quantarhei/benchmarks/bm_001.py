# -*- coding: utf-8 -*-

import quantarhei as qr


def main():
    
    with qr.energy_units("1/cm"):
        mol1 = qr.Molecule([0.0, 12000.0])
        mol2 = qr.Molecule([0.0, 12100.0])
        mol3 = qr.Molecule([0.0, 12100.0])
        
        agg = qr.Aggregate([mol1, mol2, mol3])
    
        m1 = qr.Mode(100)
        mol1.add_Mode(m1)
        
        m2 = qr.Mode(100)
        mol2.add_Mode(m2)

        m3 = qr.Mode(100)
        mol3.add_Mode(m3)
        
    agg.build(mult=1)
    
    print(agg.Ntot)




