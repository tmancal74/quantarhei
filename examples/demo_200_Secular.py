# -*- coding: utf-8 -*-
"""

    Demonstration of the secularization of a relaxation tensor



"""
import numpy
import quantarhei as qr
import quantarhei.models as models

N = 3

# Arbitrary density matrix
rho = qr.ReducedDensityMatrix(dim = N)
rho.data[0,0] = 0.1
rho.data[1,1] = 0.8
rho.data[2,2] = 0.1
rho.data[1,2] = 0.3
rho.data[2,1] = 0.3

print(rho)

# Rate matrix

KK = qr.qm.RateMatrix(dim=N)

KK.set_rate((1,2), 1.0/100.0)
KK.set_rate((0,1), 1.0/200.0)

print("\nRate matrix")
print("===========")
print(KK.data)

###############################################################################
#
# This is how Rate matrix is applied to the density matrix
#
rho_out = numpy.diagflat(numpy.einsum("ij,jj", KK.data, rho.data))
#
###############################################################################
#
# It is equivalent to the code below
#
rho_check = numpy.zeros((N,N), dtype=qr.COMPLEX)
for ii in range(N):
    for kk in range(N):
        rho_check[ii,ii] += KK.data[ii,kk]*rho.data[kk,kk]

# Checking results
        
print("\nResult of einsum:")  
print(rho_out)
print("Result of a check:")
print(rho_check)
print("Max difference:", numpy.amax(numpy.abs(rho_out-rho_check)))

# Pure dephasing from rates

GG = numpy.zeros((N,N), dtype=qr.REAL)

for ii in range(N):
    for jj in range(N):
        if ii != jj:
            GG[ii,jj] = -(KK.data[ii,ii]+KK.data[jj,jj])/2.0

print("\nDephasing matrix")
print("================")
print(GG)

###############################################################################
#
#  This is how dephasing matrix is applied to the density matrix
#
rho_out = -GG*rho.data
#
###############################################################################
#
# It is equivalent to the code below
#
rho_check = numpy.zeros((N,N), dtype=qr.COMPLEX)
for ii in range(N):
    for jj in range(N):
        rho_check[ii,jj] = -GG[ii,jj]*rho.data[ii,jj]

# checking results

print("\nResult of matrix operations:")  
print(rho_out)
print("Result of a check:")
print(rho_check)
print("Max difference:", numpy.amax(numpy.abs(rho_out-rho_check)))


RR = numpy.zeros((N,N,N,N), dtype=qr.COMPLEX)
RR[1,1,2,2] = 2.0
RR[0,0,1,1] = 1.0
RR[1,2,1,2] = -0.5
RR[1,1,1,1] = -1.0

rr = numpy.einsum("iijj->ij", RR)

print(rr)
print(rr[1,2])
print(rr[0,1])
print(rr[1,1])

gg = numpy.einsum("ijij->ij", RR)
for ii in range(gg.shape[0]):
    gg[ii,ii] = 0.0
    
print(gg)
print(gg[1,2])


R = qr.qm.RelaxationTensor()

R.data = RR

R.secularize(legacy=False)
print(R.secular_KK)
print(R.secular_GG)

mg = models.ModelGenerator()

timeaxis = qr.TimeAxis(0.0, 100, 1.0)
agg = mg.get_Aggregate_with_environment(name="trimer-1_env",
                                        timeaxis=timeaxis)
agg.build()

(Rfld, ham) = agg.get_RelaxationTensor(timeaxis,"standard_Redfield")
ham = agg.get_Hamiltonian()
sbi = agg.get_SystemBathInteraction()
#
#ham.protect_basis()
#with qr.eigenbasis_of(ham):
#    Rfld = qr.qm.RedfieldRelaxationTensor(ham, sbi, as_operators=False)
#    
#ham.unprotect_basis()
with qr.eigenbasis_of(ham):    
    Rfld.secularize(legacy=False)
print("Dephasing from Redfield tensor")
print(Rfld.secular_GG)
print("Rate matrix from Redfield tensor")
print(numpy.real(Rfld.secular_KK))

Rfld_M = qr.qm.RedfieldRateMatrix(ham, sbi)
with qr.eigenbasis_of(ham):
    print(Rfld_M.data)