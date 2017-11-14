import numpy
cimport numpy 

def loopit(numpy.ndarray Km, numpy.ndarray Kd, numpy.ndarray Lm, 
           numpy.ndarray Ld, int Na, numpy.ndarray RR, int m):

    cdef int a, b, c, d
    cdef numpy.ndarray KdLm = numpy.zeros((Na, Na), dtype=numpy.float)
    cdef numpy.ndarray LdKm = numpy.zeros((Na, Na), dtype=numpy.float)
    
    KdLm = numpy.dot(Kd,Lm[m,:,:])
    LdKm = numpy.dot(Ld[m,:,:],Km[m,:,:])

    for a in range(Na):
        for b in range(Na):
            for c in range(Na):
                for d in range(Na):
                    
                    RR[a,b,c,d] += (Km[m,a,c]*Ld[m,d,b] 
                                    + Lm[m,a,c]*Kd[d,b])
                    if b == d:
                        RR[a,b,c,d] -= KdLm[a,c] 
                    if a == c:
                        RR[a,b,c,d] -= LdKm[d,b]
            