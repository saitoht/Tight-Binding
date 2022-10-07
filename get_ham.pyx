cimport cython
import numpy as np
cimport numpy as cnp

@cython.boundscheck(False)
@cython.wraparound(False)
def get_hamkij(int nwf, tr, rv, k): 
    def ekr(k, r):
        cdef double complex kdotr, ekr
        kdotr = np.dot(np.array(k), np.array(r))
        ekr = np.cos(kdotr) + 1.j * np.sin(kdotr)
        return ekr

    cdef long n, m, l, rlen
    cdef double complex Hk_comp, phase
    cdef cnp.ndarray[cnp.complex128_t,ndim=2] Hk
    rlen = nwf * nwf
    rv = np.array(rv)
    Hk = np.zeros((nwf,nwf), dtype=np.complex128)
    for n in range(nwf):
        for m in range(nwf):
            Hk_comp = 0.+0.j
            for l in range(rlen):
                phase = ekr(k,rv[:,l])
                Hk_comp = Hk_comp + tr[l,m,n] * phase
            Hk[m,n] = Hk_comp
    return Hk
  
