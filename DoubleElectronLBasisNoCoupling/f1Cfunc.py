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
from indices import indgen
import sys
#from pymatsolver import Pardiso

R = [5.0e4, 5.0e4, 1.0e5] # [5.0e8, 5.0e8, 1.0e9]
g = [2.0087, 2.0057, 2.0021] #g factors
A = [6, 6, 36] #Gauss
B0 = 3300 #3300 #Gauss
g0 = 2.0021 #free e g factor

#diffusion tilt (all angles in degrees)
ald = 0#-5 
bed = 0#10
gad = 0#-15

#magnetic tilt
alm = 0#-20
bem = 0#25
gam = 0#-30

#dipolar tilt !!!NEW!!!
aldip = 0#-20
bedip = 0#25
gadip = 0#-30

#director tilt
psi = 0

#nuclear spin
In = 1 #nitroxide ideally
I2 = 2 * In
ipnmx = 2 

#high field approx must be valid
if (B0 < 10 * max( max(A), (max(g)-min(g))*B0/g0) ):
    raise "High field approximation violated it seems"

#generate the matrix indices, diag and off-diag spaces   
print("\n")
ndimo, ndimd, ind_offdiag, ind_diag, Lstarts_offdiag, Lstarts_diag, Llist = indgen(75, 51, 25, 2, 2, g, A, I2, [ald,bed,gad], [alm,bem,gam], [aldip, bedip, gadip], psi) #(75, 51, 25, 2, 2) #(30, 21, 15, 10, 2) #(8, 7, 4, 3, 2)
print(Llist)
print(Lstarts_diag)
np.savetxt('ind_offdiag.txt',ind_offdiag,fmt='%d')
np.savetxt('ind_diag.txt',ind_diag,fmt='%d')
#from dlkmo import *
stt = time.time()
#parms being copied from the first cell

iden = coo_matrix((np.ones(ndimo),(np.arange(ndimo),np.arange(ndimo))))
stvx = np.array([[1.],[1.],[1.]])/math.sqrt(3);stvx=np.concatenate((stvx,np.zeros((ndimo-3,1))),axis=0)

matzi = SpinHamiltonianSuperoperatorDiag(ind_diag, B0, psi, In, g, A, ald,bed,gad, alm,bem,gam, Llist, Lstarts_diag)
matzr = DiffTensorNoPotential(ind_diag, R, g)
matz = matzr + 1.0j * matzi
io.mmwrite('matz',matz.tocoo())

matxi = SpinHamiltonianSuperoperatorOffDiag(ind_offdiag, B0, psi, In, g, A, ald,bed,gad, alm,bem,gam, Llist, Lstarts_offdiag) - B0*iden
matxr = DiffTensorNoPotential(ind_offdiag, R, g)
matx = matxr + 1.0j * matxi 
io.mmwrite('matx',matx.tocoo())

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

