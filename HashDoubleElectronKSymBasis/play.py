import sys
from scipy.sparse.linalg import LinearOperator,gmres
from scipy import *
import numpy as np
def mv(v):
    return array([ 2*v[0], 3*v[1]])
A = LinearOperator( (2,2), matvec=mv )
#A.matvec( ones(2) )
#A * ones(2)
print(gmres(A,np.ones(2)))
sys.exit()
from scipy import matrix
from scipy.sparse.linalg import gmres
M = matrix([[1,2,3],[4,5,6],[7,8,25]]).tocsr()
print(gmres(aslinearoperator(M),np.ones((3,1))))
