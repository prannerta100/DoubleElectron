import numpy as np
import csv
x = list(csv.reader(open("matdip_Plus1.txt", "r"), delimiter=","))
#x = list(csv.reader(open("matzi.txt", "r"), delimiter=","))
#x = list(csv.reader(open("matxi.txt", "r"), delimiter=","))
x = np.array(x).astype(float)
I = x[:,0].astype(int)
J = x[:,1].astype(int)
E = x[:,2]+1.0j*x[:,3]
import scipy.sparse
m = scipy.sparse.coo_matrix((E,(I,J)))


from scipy.io import mmread
f=mmread('matdip_Plus1.mtx')
#f=mmread('matz.mtx')
#f=mmread('matx.mtx')

#error between f and m, should be 0
err=f-m
print('f norm',np.max(np.abs(f)),' f.shape',f.shape)
print('m norm',np.max(np.abs(m)),' m.shape',m.shape)
print('err norm',np.max(np.abs(err)))
