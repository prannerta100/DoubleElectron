import numpy as np
import scipy.sparse
import csv
def conv2coo(filename,dim, delimiter=','):
    #x = list(csv.reader(open(filename, "r"), delimiter=delimiter))
    #x = np.array(x).astype(float)
    x = np.loadtxt(filename, delimiter=delimiter)
    I = x[:,0].astype(int)
    J = x[:,1].astype(int)
    if np.max(I) > dim-1: #Fortran style array
        I -= 1
    if np.max(J) > dim-1: #Fortran style array
        J -= 1
    E = x[:,2]
    if x.shape[1] > 3: #in case imaginary elements are around
        E = E.astype(complex)
        E += 1.0j*x[:,3] #add the last column if it exists
        
    return scipy.sparse.coo_matrix((E,(I,J)),shape=(dim,dim))
