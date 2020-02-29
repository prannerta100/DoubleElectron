import numpy as np
from scipy.sparse import coo_matrix
with open('matlab_i.mtxx') as file:
    array2d = [[digit for digit in line.split()] for line in file]
I = [0] * len(array2d)
J = [0] * len(array2d)
V = [0] * len(array2d)
for i in range(len(array2d)):
    I[i], J[i], V[i] = int(array2d[i][0])-1, int(array2d[i][1])-1, float(array2d[i][2])
nlsl = -coo_matrix((V,(I,J))) #because the code does \Gamma "-" i L

with open('matrx_imag.mtxx0') as file:
    array2d = [[digit for digit in line.split()] for line in file]
I = [0] * len(array2d)
J = [0] * len(array2d)
V = [0] * len(array2d)
for i in range(len(array2d)):
    I[i], J[i], V[i] = int(array2d[i][0])-1, int(array2d[i][1])-1, float(array2d[i][2])
nlspmc = -coo_matrix((V,(I,J))) #because the code does \Gamma "-" i L

print('max diff between nlsl and nlspmc',np.max(np.abs(nlsl-nlspmc)))
