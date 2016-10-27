# -*- coding: utf-8 -*-
from quantarhei.core.matrixdata import MatrixData


file = TemporaryFile()

m = MatrixData(data=[0.0, 1.0, 2.0])
print(m.data)

m.save_data("/Users/tomas/test.dat")

m2 = MatrixData()
m2.load_data("/Users/tomas/test.dat")
print(m2.data)

m.save_data("/Users/tomas/test.txt")

m3 = MatrixData()
m3.load_data("/Users/tomas/test.dat")
print(m3.data)

m.save_data("/Users/tomas/test.npy")

m3 = MatrixData()
m3.load_data("/Users/tomas/test.npy")
print(m3.data)

m.save_data("/Users/tomas/test.npz")

m4 = MatrixData()
m4.load_data("/Users/tomas/test.npy")
print(m4.data)

