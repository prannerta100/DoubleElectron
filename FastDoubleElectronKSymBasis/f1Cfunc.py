#check nlspmc ndimo, ndimd
#NUCLEAR ZEEMAN IS 0!!!
#words of caution: we are dealing with very big matrices here, so memory allocation must be parsimonious
#no dense matrices!!!! Also, it seems adding an (n,1) and an (n,) matrix requires the generation of a (n,n) matrix. Here n~100000!
#usual caution: Python is oop, so use b = a.deepcopy() rather than b = a
from __future__ import division
import numpy as np
import math
from scipy.sparse import coo_matrix , csr_matrix
from scipy.sparse.linalg import gmres, expm_multiply, LinearOperator
from stvec_gen import stvec_gen
#from pymatsolver import Pardiso
from scipy import io
import scipy
#from myModule import w3jmatlabC, dlkmC, Ax_aOffDiagDiffC, Ax_aOffDiagSumC, Ax_gOffDiagDiffC, Ax_gOffDiagSumC, Ax_aDiagDiffC
#from SpinHamiltonianSuperoperator import SpinHamiltonianSuperoperatorDiag, SpinHamiltonianSuperoperatorOffDiag
#from DiffusionSuperoperator import DiffTensorNoPotential
#from DipolarSuperoperator import SecularDipolarTensor
import time
from common_defs import *
from indices import indgen
import os,sys
from expokit import expv
from convert_mat import conv2coo
#from pymatsolver import Pardiso

#R = [5.0e4, 5.0e4, 1.0e5] 
R = [5.0e8, 5.0e8, 1.0e9]
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
psi = 0.00001

#nuclear spin
In = 1 #nitroxide ideally
I2 = 2 * In
ipnmx = 2 

#high field approx must be valid
if (B0 < 10 * max( max(A), (max(g)-min(g))*B0/g0) ):
    raise "High field approximation violated it seems"

#dipolar constant
D = 1 #need to change

#generate the matrix indices, diag and off-diag spaces   
print("\n")
ndimo, ndimd, ndimPlus1, ndimPlus2, ndimMinus1, ind_offdiag, ind_diag, ind_Plus1, ind_Plus2, ind_Minus1, \
Lstarts_offdiag, Lstarts_diag, Lstarts_Plus1, Lstarts_Plus2, Lstarts_Minus1, Llist = \
indgen(30, 21, 15, 10, 2, g, A, I2, [ald,bed,gad], [alm,bem,gam], [aldip, bedip, gadip], psi) #(75, 51, 25, 2, 2) #(30, 21, 15, 10, 2) #(8, 7, 4, 3, 2)
#indgen(8, 7, 4, 3, 2, g, A, I2, [ald,bed,gad], [alm,bem,gam], [aldip, bedip, gadip], psi) #(75, 51, 25, 2, 2) #(30, 21, 15, 10, 2) #(8, 7, 4, 3, 2)
print(Llist)
print(Lstarts_diag)
numL=len(Llist)
stv = stvec_gen(ind_Plus1)
np.savetxt('ind_offdiag.txt',ind_offdiag,fmt='%d',delimiter=',')
np.savetxt('ind_diag.txt',ind_diag,fmt='%d',delimiter=',')
np.savetxt('stv.txt',stv)
np.savetxt('ind_Plus1.txt',ind_Plus1,fmt='%d',delimiter=',')
np.savetxt('ind_Plus2.txt',ind_Plus2,fmt='%d',delimiter=',')
np.savetxt('ind_Minus1.txt',ind_Minus1,fmt='%d',delimiter=',')
np.savetxt('Llist.txt',Llist,fmt='%d')
np.savetxt('Lstarts_offdiag.txt',Lstarts_offdiag,fmt='%d')
np.savetxt('Lstarts_diag.txt',Lstarts_diag,fmt='%d')
np.savetxt('Lstarts_Plus1.txt',Lstarts_Plus1,fmt='%d',delimiter=',')
np.savetxt('Lstarts_Plus2.txt',Lstarts_Plus2,fmt='%d',delimiter=',')
np.savetxt('Lstarts_Minus1.txt',Lstarts_Minus1,fmt='%d',delimiter=',')
os.system('gcc -o gen_mat gen_mat.c -lm ')
os.system('./gen_mat '+str(B0)+' '+str(psi)+' '+str(In)+' '+str(g[0])+' '+str(g[1])+' '+str(g[2])+' '+str(A[0])+' '+str(A[1])+\
          ' '+str(A[2])+' '+str(R[0])+' '+str(R[1])+' '+str(R[2])+' '+str(ald)+' '+str(bed)+' '+str(gad)+' '+str(alm)+\
           ' '+str(bem)+' '+str(gam)+' '+str(aldip)+' '+str(bedip)+' '+str(gadip)+' '+str(D)+' '+str(ndimo)+' '+str(ndimd)+\
           ' '+str(ndimPlus1)+' '+str(ndimPlus2)+' '+str(ndimMinus1)+' '+str(numL))
sys.exit()
stt = time.time()
matzi = conv2coo('matzi.txt')
matzr = DiffTensorNoPotential(ind_diag, R, g)
matz = matzr + 1.0j * matzi

matxi = conv2coo('matxi.txt');#B0 * iden already subtracted
matxr = DiffTensorNoPotential(ind_offdiag, R, g)
matx = matxr + 1.0j * matxi

matdipPlus1 = 1.0j * conv2coo('matdip_Plus1.txt')

def matvec_Plus1(v):
    #v expected to be in the pS1+pS2 = 1 coherence subspace 
    #reshape v in 2 ways, one for electron 1, other for electron 2
    #indx1 = np.logical_and(np.logical_and(p1arr[:,8]==2,p1arr[:,9]==0),np.logical_and(p1arr[:,10]==1,p1arr[:,11]==0))  
    v1 = np.reshape(v,(-1,36)).astype(complex) #need to have 36 columns
    offdiag_cols = [0,1,4,5,8,9,12,13,16,17,20,21,24,25,28,29,32,33] 
    diag_cols = [2,3,6,7,10,11,14,15,18,19,22,23,26,27,30,31,34,35]
    #print(matx.shape,v1.shape)
    ldh = int(len(diag_cols)/2)
    v1[:,offdiag_cols] = matx@v1[:,offdiag_cols] #off-diag cols
    #diag columns now
    tmp = np.reshape(v1[:,diag_cols],(-1,ldh,2))
    tmp = matz @ np.reshape(np.transpose(tmp,(0,2,1)),(-1,ldh)) 
    v1[:,diag_cols] = np.reshape(np.transpose(np.reshape(tmp,(-1,2,ldh)),(0,2,1)),(-1,2*ldh))
    #reshape back
    v1 = np.reshape(v1,(-1,1)) #back to the stretched out shape, Ilya style multiplication
    
    #v2 tricky, let us think this through
    #ind_Plus1 has the following hierarchy: LKMjK,(pI,qI,pS,qS)_A,(pI,qI,pS,qS)_B
    #the safe thing to do here is reshape(-1,9,9,4), then swap indices 1 and 2 among 0,1,2,3
    v2 = np.transpose(np.reshape(v,(-1,9,9,4)),(0,2,1,3))
    v2 = np.reshape(v2,(-1,36)).astype(complex)
    diag_cols = [0,1,4,5,8,9,12,13,16,17,20,21,24,25,28,29,32,33] #offdiag becomes diag and vice-versa, because 1+0=1 and 0+1=1
    offdiag_cols = [2,3,6,7,10,11,14,15,18,19,22,23,26,27,30,31,34,35]
    v2[:,offdiag_cols] = matx @ v2[:,offdiag_cols] #electron 2
    #diag columns now
    tmp = np.reshape(v2[:,diag_cols],(-1,ldh,2))
    tmp = matz @ np.reshape(np.transpose(tmp,(0,2,1)),(-1,ldh)) 
    v2[:,diag_cols] = np.reshape(np.transpose(np.reshape(tmp,(-1,2,ldh)),(0,2,1)),(-1,2*ldh))
    #reshape back
    v2 = np.reshape(v2,(-1,1)) #back to the stretched out shape, Ilya style multiplication

    #vdip, ignore for now
    vdip = matdipPlus1 @ np.reshape(v,(-1,1))
    #order 'C', a = np.arange(12) goes into 
    #array([[ 0,  1,  2],
    #       [ 3,  4,  5],
    #       [ 6,  7,  8],
    #       [ 9, 10, 11]])
    #print(v1.shape,v2.shape,vdip.shape)
    
    return v1+v2+vdip
cfact = 17
t1 = 0.004
mat1 = LinearOperator((stv.shape[0],stv.shape[0]),matvec=matvec_Plus1)
anorm = scipy.sparse.linalg.norm(matx,np.inf) + scipy.sparse.linalg.norm(matz,np.inf)
print('anorm=',anorm)
tt = expv(-cfact * t1 ,mat1, stv, anorm)
#tt = -cfact * t1 * mat1*stv
print('maxabs=',np.max(np.abs(tt[0])))
print(tt)


#tt = matvec_Plus1(stv)
#print(tt.shape,stv.shape)
#print(np.max(np.abs(tt-3*stv)))

sys.exit()

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

