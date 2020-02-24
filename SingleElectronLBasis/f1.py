#check nlspmc ndimo, ndimd
#NUCLEAR ZEEMAN IS 0!!!
from __future__ import division
import numpy as np
from sympy.physics.wigner import *
from sympy import N
def par(x):
    return 1 - 2 * (x % 2) 

g = [2.0087, 2.0057, 2.0021] #g factors
A = [6, 6, 36] #Gauss
B0 = 34050 #Gauss
g0 = 2.0021 #free e g factor

axialA = (A[0] == A[1]) #check if g axial
axialg = (g[0] == g[1]) #check if A axial

rndoff = 1e-15
#diffusion tilt (all angles in degrees)
ald = -5 
bed = 10
gad = -15
itd = abs(ald) > rndoff or abs(bed) > rndoff or abs(gad) > rndoff

#magnetic tilt
alm = -20
bem = 25
gam = -30
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

#generate the matrix    

def matogen(Lemax, Lomax, Kmax, Mmax, ipnmx):
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
ndimo, ndimd, ind_offdiag, ind_diag = matogen(8, 7, 4, 3, 2)




import scipy.special as s
import math
from rotation_matrices import dlkm
import scipy.sparse
#from dlkmo import *
#important: when I define a function within a function, I can use the variables passed to the first function
def SpinHamiltonianSuperoperatorDiag(ind_arr, B0, psi, I, g, A, ald, bed, gad, alm, bem, gam):
    np.seterr(invalid='raise')
    hbar = 1.05443e-27
    betae = 9.2731e-21
    g0 = sum(g)/3
    cfac = g0 * betae / hbar
    print(cfac)
    
    ndim = len(ind_arr)
    if len(ind_arr[0]) != 9:
        raise NameError("Problem with indices!")
    #mat = np.zeros((ndim,ndim),dtype=complex)
    mat = scipy.sparse.lil_matrix((ndim, ndim), dtype=complex)
    print(mat.shape)
    print(mat[1,2])
    
    def Wig3j(j1,j2,j3, m1,m2,m3):
        def parity(x):
            return 1-2*(x%2)
        #pass lists, not numpy arrays
        j123 = [j1,j2,j3]
        m123 = [m1,m2,m3]
        # Input error checking
        if any([x<0 for x in j123]):
            raise ValueError( 'The j must be non-negative' )
        elif any([abs(2*x-round(2*x)) > rndoff for x in j123+m123] ):
            raise ValueError( 'All arguments must be integers or half-integers' )
        #elif any( rem( (j123 - m123), 1 ) ):
        #    error( 'j123 and m123 do not match' );
        # Selection rules, j3 out of interval, non-conserving angular momentum, m is larger than j
        if j3 > j1 + j2 or j3 < abs(j1 - j2) or m1 + m2 + m3 != 0 or any([abs(m123[i]) > j123[i] for i in range(3)]): 
            return 0
        # Simple common case: #m1 = m2 = m3 = 0 & j1 + j2 + j3 is odd
        if m123 == [0,0,0] and sum(j123) % 2 == 1: 
            return 0

        # Evaluation
        t1 = j2 - m1 - j3
        t2 = j1 + m2 - j3
        t3 = j1 + j2 - j3
        t4 = j1 - m1
        t5 = j2 + m2
        tmin = max( 0,  max( t1, t2 ) )
        tmax = min( t3, min( t4, t5 ) );
        t_arr = range(tmin,tmax+1)
        arr = [j1+j2+j3+1+1, j1+j2-j3+1, j1-j2+j3+1, -j1+j2+j3+1, j1+m1+1, j1-m1+1, j2+m2+1, j2-m2+1, j3+m3+1, j3-m3+1]
        tmp = np.matmul(s.gammaln(arr), [-1,1,1,1,1,1,1,1,1,1]) * 0.5
        w = sum([ parity(t) * math.exp(-sum(s.gammaln([t+1, t-t1+1, t-t2+1, t3-t+1, t4-t+1, t5-t+1]))) for t in t_arr]) * math.exp(tmp) * parity(j1-j2-m3)
        # Warnings
        if math.isnan(w):
            print('Wigner3J is NaN!' )
        elif math.isinf(w):
            print('Wigner3J is Inf!' )
        return w

    def NormFacL(l1,l2):
        return np.sqrt((2*l1+1)*(2*l2+1))
    def NormFacK(k1,k2):
        return 1/np.sqrt(((k1==0)+1)*((k2==0)+1))
    
    f2_aa = {0:np.sqrt(2/3)*(A[2]-0.5*(A[0]+A[1])), -1:0, 1:0, -2:0.5*(A[0]-A[1]), 2:0.5*(A[0]-A[1])}
    f2_gg = {0:np.sqrt(2/3)*(g[2]-0.5*(g[0]+g[1])), -1:0, 1:0, -2:(g[0]-g[1])/2, 2:(g[0]-g[1])/2}
    d2diff = {m1:{m2:dlkm(2,m1,m2,ald,bed,gad) for m2 in range(-2,3)} for m1 in range(-2,3)}
    d2mag = {m1:{m2:dlkm(2,m1,m2,alm,bem,gam) for m2 in range(-2,3)} for m1 in range(-2,3)}
    d2diffmag = {m1:{m2:sum([d2diff[m1][mtmp] * d2mag[mtmp][m2] for mtmp in range(-2,3)]) for m2 in range(-2,3)} for m1 in range(-2,3)}
    F_aD_Conj0 = np.conj(-sum(A)/np.sqrt(3))
    F_gD_Conj0 = np.conj(-sum(g)/np.sqrt(3))
    F_aD_Conj2 = {m:sum([d2diffmag[m][mtmp]* f2_aa[mtmp].conjugate() for mtmp in range(-2,3)]) for m in range(-2,3)}
    F_gD_Conj2 = {m:sum([d2diff[m][mtmp] * f2_gg[mtmp].conjugate() for mtmp in range(-2,3)]) for m in range(-2,3)}
    #args not defined, don't want to pass too much
    #def G_a(l,jK1,jK2,K):
    #    if abs(K) > l:
    #        return 0
    #    if l==0:
    #        x = F_aD_Conj0
    #        re, im = np.real(x), -np.imag(x) #conjugate, so need to bar the number
    #    if l==2:
    #        x = F_aD_Conj2[K]
    #        re, im = np.real(x), -np.imag(x)
    #     (jK1==jK2) * re + (jK1!=jK2) * jK1 * im
    
    #def G_g(l,jK1,jK2,K):
    #    x = F_gD_Conj(l,K)
    #    re, im = np.real(x), -np.imag(x) #conjugate, so need to bar the number
    #    return (jK1==jK2) * re + (jK1!=jK2) * jK1 * im
    
    def R_a(l,jK1,jK2,L1,L2,K1,K2):
        x1 = x2 = 0 
        if l == 0:
            #x1:G_a(l,jK1,jK2,K1-K2), x2:G_a(l,jK1,jK2,K1+K2)
            if K1-K2 == 0:
                x1 = (jK1==jK2) * np.real(F_aD_Conj0) + (jK1!=jK2) * jK1 * -np.imag(F_aD_Conj0)
            if K1+K2 == 0:
                x2 = (jK1==jK2) * np.real(F_aD_Conj0) + (jK1!=jK2) * jK1 * -np.imag(F_aD_Conj0)
        if l == 2:
            #x1:G_a(l,jK1,jK2,K1-K2), x2:G_a(l,jK1,jK2,K1+K2)
            if abs(K1-K2) <=2:
                x1 = (jK1==jK2) * np.real(F_aD_Conj2[K1-K2]) + (jK1!=jK2) * jK1 * -np.imag(F_aD_Conj2[K1-K2])
            if abs(K1+K2) <=2:
                x2 = (jK1==jK2) * np.real(F_aD_Conj2[K1+K2]) + (jK1!=jK2) * jK1 * -np.imag(F_aD_Conj2[K1+K2]) 
        return (Wig3j(L1,l,L2,K1,K2-K1,-K2)) * x1 + jK2 * par(L2+K2) * (Wig3j(L1,l,L2,K1,-K2-K1,K2)) * x2
    
    def R_g(l,jK1,jK2,L1,L2,K1,K2):
        x1 = x2 = 0 
        if l == 0:
            #x1:G_g(l,jK1,jK2,K1-K2), x2:G_g(l,jK1,jK2,K1+K2)
            if K1-K2 == 0:
                x1 = (jK1==jK2) * np.real(F_gD_Conj0) + (jK1!=jK2) * jK1 * -np.imag(F_gD_Conj0)
            if K1+K2 == 0:
                x2 = (jK1==jK2) * np.real(F_gD_Conj0) + (jK1!=jK2) * jK1 * -np.imag(F_gD_Conj0)
        if l == 2:
            #x1:G_g(l,jK1,jK2,K1-K2), x2:G_g(l,jK1,jK2,K1+K2)
            if abs(K1-K2) <=2:
                x1 = (jK1==jK2) * np.real(F_gD_Conj2[K1-K2]) + (jK1!=jK2) * jK1 * -np.imag(F_gD_Conj2[K1-K2])
            if abs(K1+K2) <=2:
                x2 = (jK1==jK2) * np.real(F_gD_Conj2[K1+K2]) + (jK1!=jK2) * jK1 * -np.imag(F_gD_Conj2[K1+K2]) 
        return (Wig3j(L1,l,L2,K1,K2-K1,-K2)) * x1 + jK2 * par(L2+K2) * (Wig3j(L1,l,L2,K1,-K2-K1,K2)) * x2
        
    
    
    
    
    def R_g(l,jK1,jK2,L1,L2,K1,K2):
        return (Wig3j(L1,l,L2,K1,K2-K1,-K2)) * G_g(l,jK1,jK2,K1-K2) +\
               jK2 * par(L2+K2) * (Wig3j(L1,l,L2,K1,-K2-K1,K2)) * G_g(l,jK1,jK2,K1+K2)
    
    #assuming that the delta(p) in the expressions is pS1-pS2 + pI1-pI2, i.e., delta(pS)+delta(qS)
    def Ax_g(l,pI1,pI2,qI1,qI2):
        #assuming that the m in Ax(l,m) already satisfies m = pI1 - pI2, so this is not a general expression for Ax(l,m)
        return 0 #(no B0 term in diag space)
    def Ax_a(l,pI1,pI2,qI1,qI2):
        dpI = pI1-pI2
        dqI = qI1-qI2
        #assuming that the m in Ax(l,m) already satisfies m = pI1 - pI2, so this is not a general expression for Ax(l,m)
        tmp = qI1 * dqI + pI1 * dpI
        K_I = np.sqrt(complex(I*(I+1) - tmp*(tmp-2)/4))
        S_A = (dpI == 0) * pI1/2 + (dpI != 0) * -dqI * K_I/np.sqrt(8)
        return (abs(dpI) == abs(dqI)) * par(dpI) * np.sqrt(2*l+1) * (Wig3j(1,1,l,0,dpI,-dpI)) * S_A
    print(ald,bed,gad,alm,bem,gam, g, A)
    
    
    #Unit tests:
    if abs(NormFacL(3,5) - np.sqrt(77)) > rndoff:
        raise ValueError("NormFacL unit test failed")
    if abs(NormFacK(0,0) - 0.5) > rndoff or abs(NormFacK(0,2) - 1/np.sqrt(2)) > rndoff or abs(NormFacK(3,0) - 1/np.sqrt(2)) > rndoff or abs(NormFacK(4,5) - 1) > rndoff:
        raise ValueError("NormFacK unit test failed")
    ald_orig,bed_orig,gad_orig,alm_orig,bem_orig,gam_orig, g_orig, A_orig = ald,bed,gad,alm,bem,gam, g, A
    ald,bed,gad,alm,bem,gam, g, A = 5,10,15, 20,25,30, [complex(0,6),7,complex(36,71)], [complex(0,6),7,complex(36,71)]
    d2diff = {m1:{m2:dlkm(2,m1,m2,ald,bed,gad) for m2 in range(-2,3)} for m1 in range(-2,3)}
    d2mag = {m1:{m2:dlkm(2,m1,m2,alm,bem,gam) for m2 in range(-2,3)} for m1 in range(-2,3)}
    d2diffmag = {m1:{m2:sum([d2diff[m1][mtmp] * d2mag[mtmp][m2] for mtmp in range(-2,3)]) for m2 in range(-2,3)} for m1 in range(-2,3)}
    f2_aa = {0:np.sqrt(2/3)*(A[2]-0.5*(A[0]+A[1])), -1:0, 1:0, -2:0.5*(A[0]-A[1]), 2:0.5*(A[0]-A[1])}
    f2_gg = {0:np.sqrt(2/3)*(g[2]-0.5*(g[0]+g[1])), -1:0, 1:0, -2:(g[0]-g[1])/2, 2:(g[0]-g[1])/2}
    F_aD_Conj0 = np.conj(-sum(A)/np.sqrt(3))
    F_gD_Conj0 = np.conj(-sum(g)/np.sqrt(3))
    F_aD_Conj2 = {m:sum([d2diffmag[m][mtmp]* f2_aa[mtmp].conjugate() for mtmp in range(-2,3)]) for m in range(-2,3)}
    F_gD_Conj2 = {m:sum([d2diff[m][mtmp] * f2_gg[mtmp].conjugate() for mtmp in range(-2,3)]) for m in range(-2,3)}

    if abs(F_gD_Conj0 - -(43-77j)/np.sqrt(3)) > rndoff or abs(F_gD_Conj2[1]-(-5.311264597718748+11.996830758823908j)) > 1e-10:
        raise ValueError("F_gD_Conj unit test failed")
    if abs(F_aD_Conj0 - -(43-77j)/np.sqrt(3)) > rndoff or abs(F_aD_Conj2[1]-(2.50766627945+37.0389445955j)) > 1e-10:
        raise ValueError("F_aD_Conj unit test failed")
    #if abs(G_a(0,1,1,0) - -43/np.sqrt(3)) > rndoff or abs(G_a(0,-1,4,0) - 77/np.sqrt(3)) > rndoff:
    #    raise ValueError("G_a unit test failed")
    #if abs(G_g(0,1,1,0) - -43/np.sqrt(3)) > rndoff or abs(G_g(0,-1,4,0) - 77/np.sqrt(3)) > rndoff:
    #    raise ValueError("G_g unit test failed")
    if abs(Ax_g(2,1,1,0,0)-0) > rndoff:
        raise ValueError("Ax_g unit test failed")
    if abs(Ax_a(2,1,0,0,1) - 3/8) > rndoff: #nuclear spin flip happening
        raise ValueError("Ax_a unit test failed")
    print("Great, all unit tests passed")
    ald,bed,gad,alm,bem,gam, g, A = ald_orig,bed_orig,gad_orig,alm_orig,bem_orig,gam_orig, g_orig, A_orig
    d2diff = {m1:{m2:dlkm(2,m1,m2,ald,bed,gad) for m2 in range(-2,3)} for m1 in range(-2,3)}
    d2mag = {m1:{m2:dlkm(2,m1,m2,alm,bem,gam) for m2 in range(-2,3)} for m1 in range(-2,3)}
    d2diffmag = {m1:{m2:sum([d2diff[m1][mtmp] * d2mag[mtmp][m2] for mtmp in range(-2,3)]) for m2 in range(-2,3)} for m1 in range(-2,3)}
    f2_aa = {0:np.sqrt(2/3)*(A[2]-0.5*(A[0]+A[1])), -1:0, 1:0, -2:0.5*(A[0]-A[1]), 2:0.5*(A[0]-A[1])}
    f2_gg = {0:np.sqrt(2/3)*(g[2]-0.5*(g[0]+g[1])), -1:0, 1:0, -2:(g[0]-g[1])/2, 2:(g[0]-g[1])/2}
    F_aD_Conj0 = np.conj(-sum(A)/np.sqrt(3))
    F_gD_Conj0 = np.conj(-sum(g)/np.sqrt(3))
    F_aD_Conj2 = {m:sum([d2diffmag[m][mtmp]* f2_aa[mtmp].conjugate() for mtmp in range(-2,3)]) for m in range(-2,3)}
    F_gD_Conj2 = {m:sum([d2diff[m][mtmp] * f2_gg[mtmp].conjugate() for mtmp in range(-2,3)]) for m in range(-2,3)}
    
    print(ndim,'x',ndim, 'matrix')
    print(ald,bed,gad,alm,bem,gam, g, A)
    
    #stt = time.time()
    #for i in range(ndim):
    #    L1,M1,K1,jM1,jK1,pI1,qI1,pS1,qS1 = ind_arr[i]
    #    for j in range(i+1):
    #        L2,M2,K2,jM2,jK2,pI2,qI2,pS2,qS2 = ind_arr[j]
    #        xcv = Wig3j(8,7,5,0,1,-1)
    #stp=time.time()
    #print(stp-stt)
    
    for i in range(ndim):
        L1,M1,K1,jM1,jK1,pI1,qI1,pS1,qS1 = ind_arr[i]
        #print(L1,M1,K1,jM1,jK1,pI1,qI1,pS1,qS1)
        #if pS1 != 0 or qS1 != 1:
        #    raise ValueError("wrong routine used, check index list again")
            
        for j in range(i+1):
            L2,M2,K2,jM2,jK2,pI2,qI2,pS2,qS2 = ind_arr[j]
            #print(L1,M1,K1,jM1,jK1,pI1,qI1,pS1,qS1)
            #if pS2 != 0 or qS2 != 1:
            #    raise ValueError("wrong routine used, check index list again")
            dpI = pI1 - pI2
            dqI = qI1 - qI2
            #neglecting nuclear Zeeman 
            #check these conditions again, I put them to speed up, may not work if initial settings change
            if abs(dpI) <=1 and abs(dpI)==abs(dqI) and abs(M2-M1) <=2:
                mat[i,j] = NormFacL(L1,L2) * NormFacK(K1,K2) * par(M1+K1) * \
                        sum([R_a(l,jK1,jK2,L1,L2,K1,K2) * (Wig3j(L1,l,L2,M1,M2-M1,-M2)) * \
                             Ax_a(l,pI1,pI2,qI1,qI2) * dlkm(l,dpI,M1-M2,0,psi,0) for l in [0,2]]) 
                mat[j,i] = mat[i,j]
    
    
    return mat
import time
stt = time.time()
#parms being copied from the first cell
matxi = SpinHamiltonianSuperoperatorDiag(ind_diag, B0, psi, In, [2.0087,2.0057,2.0021], [6,6,36], ald,bed,gad, alm,bem,gam)
#matxr = 
stp = time.time()
print(stp-stt)
with open('matrx_imag.mtxz') as file:
    array2d = [[digit for digit in line.split()] for line in file]
I = [0] * len(array2d)
J = [0] * len(array2d)
V = [0] * len(array2d)
for i in range(len(array2d)):
    I[i], J[i], V[i] = int(array2d[i][0])-1, int(array2d[i][1])-1, float(array2d[i][2])
a = -scipy.sparse.coo_matrix((V,(I,J))) #because the code does \Gamma "-" i L
print(a.shape)
#print max deviation b/w nlspmc and this routine
#print(np.max(np.abs(a-matxi)))
