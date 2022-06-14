# -*- coding: utf-8 -*-
import numpy
import matplotlib.pylab as plt

from quantarhei import TimeAxis
from quantarhei.qm import RateMatrix
from quantarhei import PopulationPropagator

tt = TimeAxis(0.0, 5000, 0.1)

N = 4 # number of sites

#
# Rate matrix (definition and initialization)
#
KK = RateMatrix(data=numpy.zeros((N,N),dtype=numpy.float64))

T01 = 25.0
T12 = 10.0
T23 = 20.0
T02 = 50.0

f01 = 0.10
f12 = 0.01
f23 = 0.1
f02 = 0.001

KK.set_rate((2,3),1.0/T23)
KK.set_rate((3,2),f23/T23)
KK.set_rate((0,1),1.0/T01)
KK.set_rate((1,0),f01/T01)
KK.set_rate((1,2),1.0/T12)
KK.set_rate((2,1),f12/T12)
KK.set_rate((0,2),1.0/T02)
KK.set_rate((2,0),f02/T02)
print("")
print("Rate matrix")
print(KK.data)


#
# initial population
#
pini = numpy.zeros(N,dtype=numpy.float64)
pini[3] = 1.0
pini[2] = 3.0
pini[1] = 2.0
pini[0] = 1.0
pini = pini/7.0


# Population propagator can be created using KK.data (matrix) or KK (RateMatrix object)
pp = PopulationPropagator(tt,KK)


#
# Exact propagated populations
#
pops_exact = pp.propagate(pini)

#
# Propagation matrix
#
Ue,(Uc0e,Uc1e,Uc2e) = pp.get_PropagationMatrix(tt, corrections=2, exact=True)
U,cors = pp.get_PropagationMatrix(tt, corrections=8, exact=False)
Uc0 = cors[0]
Uc1 = cors[1]
Uc2 = cors[2]
Uc3 = cors[3]
Uc4 = cors[4]
Uc5 = cors[5]
Uc6 = cors[6]
Uc7 = cors[7]
Uc8 = cors[8]

S2 = Uc0+Uc1+Uc2
S2e = Uc0e+Uc1e+Uc2e
S3 = S2+Uc3 #+Uc4+Uc5
S4 = S3+Uc4
S8 = S4+Uc5+Uc6+Uc7+Uc8
D2 = U - S2

k = 0
pops_U = numpy.zeros((tt.length,N),dtype=numpy.float64)
pops_Uc0 = numpy.zeros((tt.length,N),dtype=numpy.float64)
pops_Uc1 = numpy.zeros((tt.length,N),dtype=numpy.float64)
pops_Uc1e = numpy.zeros((tt.length,N),dtype=numpy.float64)
pops_Uc2 = numpy.zeros((tt.length,N),dtype=numpy.float64)
pops_Uc2e = numpy.zeros((tt.length,N),dtype=numpy.float64)
pops_Uc3 = numpy.zeros((tt.length,N),dtype=numpy.float64)
pops_Uc4 = numpy.zeros((tt.length,N),dtype=numpy.float64)
pops_Uc8 = numpy.zeros((tt.length,N),dtype=numpy.float64)
pops_D = numpy.zeros((tt.length,N),dtype=numpy.float64)

for t in tt.data:
    pops_U[k,:] = numpy.dot(U[:,:,k],pini)
    pops_Uc0[k,:] = numpy.dot(Uc0[:,:,k],pini)
    pops_Uc1[k,:] = numpy.dot(Uc0[:,:,k]+Uc1[:,:,k],pini)
    pops_Uc1e[k,:] = numpy.dot(Uc0e[:,:,k]+Uc1e[:,:,k],pini)
    pops_Uc2[k,:] = numpy.dot(S2[:,:,k],pini)
    pops_Uc2e[k,:] = numpy.dot(S2e[:,:,k],pini)
    pops_Uc3[k,:] = numpy.dot(S3[:,:,k],pini)
    pops_Uc4[k,:] = numpy.dot(S4[:,:,k],pini)
    pops_Uc8[k,:] = numpy.dot(S8[:,:,k],pini)
    pops_D[k,:] = numpy.dot(D2[:,:,k],pini)
    k += 1


#
# Plotting populations
#
#plt.plot(tt.data,pops_exact[:,2])
#plt.plot(tt.data,pops_exact[:,1])
#plt.plot(tt.data,pops_exact[:,0])

#plt.show()

plt.plot(tt.data,pops_U[:,3],"-k")
plt.plot(tt.data,pops_U[:,2],"-b")
plt.plot(tt.data,pops_U[:,1],"-r")
plt.plot(tt.data,pops_U[:,0],"-g")
#plt.plot(tt.data,pops_Uc0[:,3],"--k")
#plt.plot(tt.data,pops_Uc0[:,2],"--b")
#plt.plot(tt.data,pops_Uc0[:,1],"--r")
#plt.plot(tt.data,pops_Uc0[:,0],"--g")
plt.plot(tt.data,pops_Uc1e[:,3],"-k")
plt.plot(tt.data,pops_Uc1e[:,2],"-b")
plt.plot(tt.data,pops_Uc1e[:,1],"-r")
plt.plot(tt.data,pops_Uc1e[:,0],"-g")
plt.plot(tt.data,pops_Uc1[:,3],"--k")
plt.plot(tt.data,pops_Uc1[:,2],"--b")
plt.plot(tt.data,pops_Uc1[:,1],"--r")
plt.plot(tt.data,pops_Uc1[:,0],"--g")
plt.plot(tt.data,pops_Uc2e[:,3],"-k")
plt.plot(tt.data,pops_Uc2e[:,2],"-b")
plt.plot(tt.data,pops_Uc2e[:,1],"-r")
plt.plot(tt.data,pops_Uc2e[:,0],"-g")
plt.plot(tt.data,pops_Uc2[:,3],"--k")
plt.plot(tt.data,pops_Uc2[:,2],"--b")
plt.plot(tt.data,pops_Uc2[:,1],"--r")
plt.plot(tt.data,pops_Uc2[:,0],"--g")
plt.plot(tt.data,pops_Uc3[:,3],"--k")
plt.plot(tt.data,pops_Uc3[:,2],"--b")
plt.plot(tt.data,pops_Uc3[:,1],"--r")
plt.plot(tt.data,pops_Uc3[:,0],"--g")
plt.plot(tt.data,pops_Uc4[:,3],"--k")
plt.plot(tt.data,pops_Uc4[:,2],"--b")
plt.plot(tt.data,pops_Uc4[:,1],"--r")
plt.plot(tt.data,pops_Uc4[:,0],"--g")
plt.plot(tt.data,pops_Uc8[:,3],"--k")
plt.plot(tt.data,pops_Uc8[:,2],"--b")
plt.plot(tt.data,pops_Uc8[:,1],"--r")
plt.plot(tt.data,pops_Uc8[:,0],"--g")
plt.show()

#
# Markovian remainder
#
#plt.plot(tt.data,pops_D[:,2])
#plt.plot(tt.data,pops_D[:,1])
#plt.plot(tt.data,pops_D[:,0])

#plt.show()

#
#  Plotting propagation matrix
#
#for i in range(N):
#    for j in range(N):
#        plt.plot(tt.data,U[i,j,:])
#plt.show()
    
#for i in range(N):
#    for j in range(N):
#        plt.plot(tt.data,Uc0[i,j,:]+Uc1[i,j,:],"--")
        
#plt.show()

