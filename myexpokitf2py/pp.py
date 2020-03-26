from expokit import expv
import scipy
from numpy import squeeze
import numpy as np
from scipy.sparse.linalg import LinearOperator
stv = scipy.ones((5,1))
print(stv.shape[0])
def mv(v):
    y = squeeze(v)
    return scipy.array([y[0]+2*y[1],3*y[0]+4*y[1],y[2],y[3],y[4]])
mat1 = LinearOperator((stv.shape[0],stv.shape[0]),matvec=mv)
#print(np.linalg.norm(mat1))
#print(expv(10,scipy.eye(4),scipy.ones((4,1))))
print(expv(1,mat1,stv,7))

