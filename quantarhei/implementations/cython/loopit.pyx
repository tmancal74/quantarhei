import numpy
cimport numpy 

REAL_t = numpy.float64_t

def loopit(numpy.ndarray Km, numpy.ndarray Kd, numpy.ndarray Lm, 
           numpy.ndarray Ld, long Na, numpy.ndarray RR, long m):

    cdef long a, b, c, d
    cdef numpy.ndarray KdLm = numpy.zeros((Na, Na), dtype=REAL_t)
    cdef numpy.ndarray LdKm = numpy.zeros((Na, Na), dtype=REAL_t)
    
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
            