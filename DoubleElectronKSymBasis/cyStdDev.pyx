from numpy cimport ndarray, float64_t
cdef extern from "std_dev.h":
    double std_dev(double *arr, size_t siz)

def cStdDev(ndarray[float64_t, ndim=1] a not None):
    return std_dev(<double*> a.data, a.size)
