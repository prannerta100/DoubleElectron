import numpy as np
import scipy.sparse
import csv
def conv2coo(filename,dim):
    x = list(csv.reader(open(filename, "r"), delimiter=","))
    x = np.array(x).astype(float)
    I = x[:,0].astype(int)
    J = x[:,1].astype(int)
    E = x[:,2]+1.0j*x[:,3]
    return scipy.sparse.coo_matrix((E,(I,J)),shape=(dim,dim))
