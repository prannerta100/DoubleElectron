#check nlspmc ndimo, ndimd
#NUCLEAR ZEEMAN IS 0!!!
from __future__ import division
import numpy as np
import math
from scipy.sparse import coo_matrix , csr_matrix
from scipy.sparse.linalg import gmres
#from pymatsolver import Pardiso
from scipy import io
#from myModule import w3jmatlabC, dlkmC, Ax_aOffDiagDiffC, Ax_aOffDiagSumC, Ax_gOffDiagDiffC, Ax_gOffDiagSumC, Ax_aDiagDiffC
from SpinHamiltonianSuperoperator import SpinHamiltonianSuperoperatorDiag, SpinHamiltonianSuperoperatorOffDiag
from DiffusionSuperoperator import DiffTensorNoPotential
import time
from common_defs import *
import sys
#from pymatsolver import Pardiso

R = [5.0e4, 5.0e4, 1.0e5] # [5.0e8, 5.0e8, 1.0e9]
g = [2.0087, 2.0057, 2.0021] #g factors
A = [6, 6, 36] #Gauss
B0 = 3300 #3300 #Gauss
g0 = 2.0021 #free e g factor

axialA = (A[0] == A[1]) #check if g axial
axialg = (g[0] == g[1]) #check if A axial

#diffusion tilt (all angles in degrees)
ald = 0#-5 
bed = 0#10
gad = 0#-15
itd = abs(ald) > rndoff or abs(bed) > rndoff or abs(gad) > rndoff

#magnetic tilt
alm = 0#-20
bem = 0#25
gam = 0#-30
itm = abs(alm) > rndoff or abs(bem) > rndoff or abs(gam) > rndoff

#director tilt
psi = 0

#
In = 1 #nitroxide ideally
ipsi0 = (abs(psi)>rndoff) #director tilt flag; multiple ort or given director tilt
kptmx = 0 #I set it, but I don't see it defined to be 0 anywhere. They define it only when c20, etc. are on. It is argmax_K c_{LK}.
delta_K = 1 + (itd == 0 and itm == 0)
delta_L = 1 + (axialA and axialg and delta_K == 2 and kptmx == 0)
jMmin = 1 #as Zeeman Nuclear is 0
jKmin = 1 - 2 * (abs(alm) > rndoff or abs(gam) > rndoff or abs(gad) > rndoff) #as per rules from MIIMFreed, I read this in the code
ipnmx = 2 
I2 = 2 * In

#high field approx must be valid
if (B0 < 10 * max( max(A), (max(g)-min(g))*B0/g0) ):
    raise "High field approximation violated it seems"

#generate the matrix indices, diag and off-diag spaces   
def indgen(Lemax, Lomax, Kmax, Mmax, ipnmx):
    ind_offdiag = []
    ind_diag = []
    pImaxval = ipnmx
    #assuming that none of Lemax, Lomax, Kmax, Mmax are 0, and even ipnmx
    if par(Lemax) == -1:
        Lemax -= 1
    if Lomax > 0 and par(Lomax) == 1:
        Lomax -= 1
    ##if(lomx.gt.lemx): L can't cross Lemax, so no worries
    if Kmax > Lemax:
        Kmax = Lemax
    if Mmax > Lemax:
        Mmax = Lemax
    if par(Kmax) == -1:
        Kmax -= 1
    if ipnmx > I2:
        ipnmx = I2
    if ipsi0 == 0 and Mmax > ipnmx:
        Mmax = ipnmx
     
    ndimo = 0
    ndimd = 0
    for L in range(0,Lemax+1,delta_L):
        if par(L) == -1 and L > Lomax:
            continue
        for jK in range(jKmin, 1+1, 2):
            for K in range(0,min(Kmax, L)+1, delta_K):
                if K == 0 and par(L) != jK:
                    continue
                for jM in range(jMmin,1+1,2):
                    for M in range(0,min(Mmax, L)+1):
                        #calculate pImin and max
                        pImax=min(I2,pImaxval)
                        if M == 0:
                            pImin=0
                        else:
                            pImin=-pImax
                        for pI in range(pImin, pImax+1):
                            if M==0 and pI==0 and jM != par(L):
                                continue
                            if ipsi0==0 and pI != M:
                                continue
                            qImax=I2 - abs(pI)
                            qImin=-qImax
                            for qI in range(qImin, qImax+1, 2):
                                ndimo += 1
                                ind_offdiag.append([L,M,K,jM,jK,pI,qI,1,0])
                                if pI == 0 and M == 0: 
                                    ndimd += 1
                                    ind_diag.append([L,M,K,jM,jK,pI,qI,0,1])
                                else:
                                    ndimd += 2
                                    ind_diag.append([L,M,K,jM,jK,pI,qI,0,1])
                                    ind_diag.append([L,-M,K,jM,jK,-pI,qI,0,1])
    return ndimo, ndimd, ind_offdiag, ind_diag
    
print("\n")
ndimo, ndimd, ind_offdiag, ind_diag = indgen(8, 7, 4, 3, 2) #(75, 51, 25, 2, 2) #(30, 21, 15, 10, 2)

np.savetxt('ind_offdiag.txt',ind_offdiag,fmt='%d')
#from dlkmo import *
stt = time.time()
#parms being copied from the first cell

iden = coo_matrix((np.ones(ndimo),(np.arange(ndimo),np.arange(ndimo))))
stvx = np.array([[1.],[1.],[1.]])/math.sqrt(3);stvx=np.concatenate((stvx,np.zeros((ndimo-3,1))),axis=0)

matzi = SpinHamiltonianSuperoperatorDiag(ind_diag, B0, psi, In, g, A, ald,bed,gad, alm,bem,gam)
matzr = DiffTensorNoPotential(ind_diag, R, g)
matz = matzr + 1.0j * matzi

matxi = SpinHamiltonianSuperoperatorOffDiag(ind_offdiag, B0, psi, In, g, A, ald,bed,gad, alm,bem,gam) - B0*iden
matxr = DiffTensorNoPotential(ind_offdiag, R, g)
matx = matxr + 1.0j * matxi 

#time required to construct the full (re+im) SLE matrix
stp = time.time()
print(stp-stt, ' seconds')

#np.dot(stvx.T,np.matmul(matx,stvx))
#print(stvx.shape)

NGRD=512
shiftr=0.2 #Gauss
data=np.zeros((NGRD,2))
ci=(0+1j)
rnge=100 #Gauss
for j in range(NGRD):
    shift=np.complex(shiftr,-rnge+j*2*rnge/(NGRD-1)) #in Gauss
    op=matx+shift*iden
    x, info = gmres(op,stvx)
    #print('x.shape',x.shape)
    data[j][0]=-rnge+j*2*rnge/(NGRD-1)
    #print(np.vdot(stvx,x).shape)
    data[j][1]=np.real(np.vdot(stvx,x))
np.savetxt("out_spec.txt",data)

#sys.exit()
print('OffDiag space matrix')
#with open('matrx_imag.mtxx0') as file:
with open('matlab_i.mtxx') as file:
    array2d = [[digit for digit in line.split()] for line in file]
I = [0] * len(array2d)
J = [0] * len(array2d)
V = [0] * len(array2d)
for i in range(len(array2d)):
    I[i], J[i], V[i] = int(array2d[i][0])-1, int(array2d[i][1])-1, float(array2d[i][2])
a = -coo_matrix((V,(I,J))) #because the code does \Gamma "-" i L
#print(a.shape)
#print max deviation b/w nlspmc and this routine
#np.savetxt('diffmat.txt',(np.abs(a-mm)>1e-10)-1+1,fmt='%d')
#np.savetxt('mat1.txt',np.abs(mm))
#a = a - 3300 * iden
io.mmwrite('matxi.txt',matxi.tocoo())
io.mmwrite('a.txt',a.tocoo())
print(np.max(np.abs(a-matxi)))
err = a-matxi
err = err.tocoo()
print(err)
print(err.shape)
io.mmwrite('diffmat.txt',err)
print('Sym_check_matxi',np.max(np.abs(matxi-csr_matrix.transpose(matxi))))
#print(np.max(np.abs(a-3300 * iden)))

print('Diag space matrix')
with open('matrx_imag.mtxz0') as file:
    array2d = [[digit for digit in line.split()] for line in file]
I = [0] * len(array2d)
J = [0] * len(array2d)
V = [0] * len(array2d)
for i in range(len(array2d)):
    I[i], J[i], V[i] = int(array2d[i][0])-1, int(array2d[i][1])-1, float(array2d[i][2])
a = -coo_matrix((V,(I,J))) #because the code does \Gamma "-" i L
#print(a.shape)
#print max deviation b/w nlspmc and this routine
#np.savetxt('diffmat.txt',(np.abs(a-mm)>1e-10)-1+1,fmt='%d')
#np.savetxt('mat1.txt',np.abs(mm))
print(np.max(np.abs(a-matzi)))
io.mmwrite('matzi.txt',matzi.tocoo())
#print(np.max(np.imag(matzi)))

