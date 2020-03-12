#important: when I define a function within a function, I can use the variables passed to the first function
import numpy as np
import math
from LSymBasisModule import w3jmatlabC, dlkmC, Ax_aC, Ax_gC 
from scipy.sparse import lil_matrix
from common_defs import *

def SpinHamiltonianSuperoperatorDiagElectronA(ind_arr, B0, psi, I, g, A, ald, bed, gad, alm, bem, gam, Llist, Lstarts_diag):
    np.seterr(invalid='raise')
    g0 = sum(g)/3
    
    ndim = len(ind_arr)
    if len(ind_arr[0]) != 12:
        raise NameError("Problem with indices!")
    mat = lil_matrix((ndim,ndim),dtype=complex) #complex because of something, not sure what 

    d2diff = {m1:{m2:dlkmC(2,m1,m2,ald,bed,gad) for m2 in range(-2,3)} for m1 in range(-2,3)}
    d2mag = {m1:{m2:dlkmC(2,m1,m2,alm,bem,gam) for m2 in range(-2,3)} for m1 in range(-2,3)}
    d2diffmag = {m1:{m2:sum([d2diff[m1][mtmp] * d2mag[mtmp][m2] for mtmp in range(-2,3)]) for m2 in range(-2,3)} for m1 in range(-2,3)}
    f2_aa = {0:math.sqrt(2/3)*(A[2]-0.5*(A[0]+A[1])), -1:0, 1:0, -2:0.5*(A[0]-A[1]), 2:0.5*(A[0]-A[1])}
    f2_gg = {0:math.sqrt(2/3)*(g[2]-0.5*(g[0]+g[1]))/g0, -1:0, 1:0, -2:(g[0]-g[1])/(2*g0), 2:(g[0]-g[1])/(2*g0)}
    F_aD_Conj0 = np.conj(-sum(A)/math.sqrt(3))
    F_gD_Conj0 = np.conj(-sum(g)/math.sqrt(3))/g0
    F_aD_Conj2 = {m:sum([d2diffmag[m][mtmp]* f2_aa[mtmp].conjugate() for mtmp in range(-2,3)]) for m in range(-2,3)}
    F_gD_Conj2 = {m:sum([d2diff[m][mtmp] * f2_gg[mtmp].conjugate() for mtmp in range(-2,3)]) for m in range(-2,3)}

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
        #if abs(w3jmatlabC(L1,l,L2,K1,K2-K1,-K2)-Wig3j(L1,l,L2,K1,K2-K1,-K2))>1e-10:
        #    print(L1,l,L2,K1,K2-K1,-K2,Wig3j(L1,l,L2,K1,K2-K1,-K2),w3jmatlabC(L1,l,L2,K1,K2-K1,-K2))
        return (w3jmatlabC(L1,l,L2,K1,K2-K1,-K2)) * x1 + jK2 * par(L2+K2) * (w3jmatlabC(L1,l,L2,K1,-K2-K1,K2)) * x2
    
#    def Ax_g(l,pI1,pI2,qI1,qI2): 
#        #assuming that the m in Ax(l,m) already satisfies m = pI1 - pI2, so this is not a general expression for Ax(l,m)
#        return 0 #(no B0 term in diag space)

    #the delta(p) in the expressions is pS1-pS2 + pI1-pI2, i.e., delta(pS)+delta(qS)
    #print(ald,bed,gad,alm,bem,gam, g, A)
    
    print(ndim,'x',ndim, 'matrix')
    #print(ald,bed,gad,alm,bem,gam, g, A)
    
    for i in range(ndim): #ndim is ndimd
        L1,M1,K1,jK1,pI1,qI1,pS1,qS1 = ind_arr[i]
        #print(L1,M1,K1,jK1,pI1,qI1,pS1,qS1)
        #if pS1 != 0 or qS1 != 1:
        #    raise ValueError("wrong routine used, check index list again")
        #left limit
        if L1-2 in Llist:
            left_j = Lstarts_diag[Llist.index(L1-2)]
        elif L1-1 in Llist:
            left_j = Lstarts_diag[Llist.index(L1-1)]
        else:
            left_j = Lstarts_diag[Llist.index(L1)]
        #right limit
        if L1+2 in Llist:
            right_j = Lstarts_diag[Llist.index(L1+2)+1]#note the starting point right after L1+2
        elif L1+1 in Llist:
            right_j = Lstarts_diag[Llist.index(L1+1)+1]#Lstarts arrays have a length 1+len(Llist), last element must be ndim
        else:
            right_j = ndim
        
        for j in range(left_j,right_j):  #range(ndim) #define band, this will really speed up things! Think about the nbnd in CW nlsl matrll.f
            #if i==7 or j==7:
            #    print(left_j,right_j)
            L2,M2,K2,jK2,pI2,qI2,pS2,qS2 = ind_arr[j]
            #print(L1,M1,K1,jM1,jK1,pI1,qI1,pS1,qS1)
            #if pS2 != 0 or qS2 != 1:
            #    raise ValueError("wrong routine used, check index list again")
            dpI = pI1 - pI2
            dqI = qI1 - qI2
            dpS = pS1 - pS2
            dqS = qS1 - qS2
            #neglecting nuclear Zeeman 
            #check these conditions again, I put them to speed up, may not work if initial settings change
            #if abs(dpI) <=1 and abs(dpI)==abs(dqI) and abs(M2-M1) <=2: #this is from the original LMSym Lee1994 code
            if abs(M2-M1) <=2:
                mat[i,j] = NormFacL(L1,L2) * NormFacK(K1,K2) * par(M1+K1) * \
                        sum([R_a(l,jK1,jK2,L1,L2,K1,K2) * w3jmatlabC(L1,l,L2,M1,M2-M1,-M2) * \
                             Ax_aC(l,pI1,pI2,qI1,qI2,pS1,pS2,qS1,qS2,I) * dlkmC(l,dpS+dpI,M1-M2,0,psi,0) for l in [0,2]]) 
                #using Eva's 1982 paper
#                mat[j,i] = mat[i,j] #forced symmetrization, shouldn't have to use it
    
    mat = mat.tocsr()
    return mat



def SpinHamiltonianSuperoperatorOffDiagElectronA(ind_arr, B0, psi, I, g, A, ald, bed, gad, alm, bem, gam, Llist, Lstarts_offdiag):
	#mind it, this is pS=1 in the 1 electron 2d ELDOR
    np.seterr(invalid='raise')
    g0 = sum(g)/3
    
    ndim = len(ind_arr)
    if len(ind_arr[0]) != 12:
        raise NameError("Problem with indices!")
    mat = lil_matrix((ndim,ndim),dtype=complex) #keeping it complex because of something, not sure what, likely the stuff in Axa_C, etc. 

    d2diff = {m1:{m2:dlkmC(2,m1,m2,ald,bed,gad) for m2 in range(-2,3)} for m1 in range(-2,3)}
    d2mag = {m1:{m2:dlkmC(2,m1,m2,alm,bem,gam) for m2 in range(-2,3)} for m1 in range(-2,3)}
    d2diffmag = {m1:{m2:sum([d2diff[m1][mtmp] * d2mag[mtmp][m2] for mtmp in range(-2,3)]) for m2 in range(-2,3)} for m1 in range(-2,3)}
    f2_aa = {0:math.sqrt(2/3)*(A[2]-0.5*(A[0]+A[1])), -1:0, 1:0, -2:0.5*(A[0]-A[1]), 2:0.5*(A[0]-A[1])}
    f2_gg = {0:math.sqrt(2/3)*(g[2]-0.5*(g[0]+g[1]))/g0, -1:0, 1:0, -2:0.5*(g[0]-g[1])/g0, 2:0.5*(g[0]-g[1])/g0}
    F_aD_Conj0 = np.conj(-sum(A)/math.sqrt(3))
    F_gD_Conj0 = np.conj(-sum(g)/math.sqrt(3))/g0
    F_aD_Conj2 = {m:sum([d2diffmag[m][mtmp]* f2_aa[mtmp].conjugate() for mtmp in range(-2,3)]) for m in range(-2,3)}
    F_gD_Conj2 = {m:sum([d2diff[m][mtmp] * f2_gg[mtmp].conjugate() for mtmp in range(-2,3)]) for m in range(-2,3)}

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
        #if abs(w3jmatlabC(L1,l,L2,K1,K2-K1,-K2)-Wig3j(L1,l,L2,K1,K2-K1,-K2))>1e-10:
        #    print(L1,l,L2,K1,K2-K1,-K2,Wig3j(L1,l,L2,K1,K2-K1,-K2),w3jmatlabC(L1,l,L2,K1,K2-K1,-K2))
        return (w3jmatlabC(L1,l,L2,K1,K2-K1,-K2)) * x1 + jK2 * par(L2+K2) * (w3jmatlabC(L1,l,L2,K1,-K2-K1,K2)) * x2
    
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
        return (w3jmatlabC(L1,l,L2,K1,K2-K1,-K2)) * x1 + jK2 * par(L2+K2) * (w3jmatlabC(L1,l,L2,K1,-K2-K1,K2)) * x2
        
    #assuming that the delta(p) in the expressions is pS1-pS2 + pI1-pI2, i.e., delta(pS)+delta(qS)
    #print(ald,bed,gad,alm,bem,gam, g, A)
    
    print(ndim,'x',ndim, 'matrix')
    #print(ald,bed,gad,alm,bem,gam, g, A)
    
    for i in range(ndim):
        L1,M1,K1,jK1,pI1,qI1,pS1,qS1 = ind_arr[i]
        #print(L1,M1,K1,jM1,jK1,pI1,qI1,pS1,qS1)
        #if pS1 != 0 or qS1 != 1:
        #    raise ValueError("wrong routine used, check index list again")
        #left limit
        if L1-2 in Llist:
            left_j = Lstarts_offdiag[Llist.index(L1-2)]
        elif L1-1 in Llist:
            left_j = Lstarts_offdiag[Llist.index(L1-1)]
        else:
            left_j = Lstarts_offdiag[Llist.index(L1)]
        #right limit
        if L1+2 in Llist:
            right_j = Lstarts_offdiag[Llist.index(L1+2)+1]#note the starting point right after L1+2
        elif L1+1 in Llist:
            right_j = Lstarts_offdiag[Llist.index(L1+1)+1]#Lstarts arrays have a length 1+len(Llist), last element must be ndim
        else:
            right_j = ndim
            
        for j in range(left_j,right_j): #range(ndim):
            L2,M2,K2,jK2,pI2,qI2,pS2,qS2 = ind_arr[j]
            #print(L1,M1,K1,jM1,jK1,pI1,qI1,pS1,qS1)
            #if pS2 != 0 or qS2 != 1:
            #    raise ValueError("wrong routine used, check index list again")
            dpI = pI1 - pI2
            dqI = qI1 - qI2
            dpS = pS1 - pS2
            dqS = qS1 - qS2

            #neglecting nuclear Zeeman 
            #check the following conditions again, I put them to speed up, may not work if initial settings change
            i#if abs(dpI) <=1 and abs(dpI)==abs(dqI) and (abs(M1-M2) <=2 or abs(M1+M2)<=2):
            if abs(M2-M1) <=2:
                mat[i,j] = NormFacL(L1,L2) * NormFacK(K1,K2) * par(M1+K1) * (\
                        sum([R_a(l,jK1,jK2,L1,L2,K1,K2) * w3jmatlabC(L1,l,L2,M1,M2-M1,-M2) * \
                             Ax_aC(l,pI1,pI2,qI1,qI2,pS1,pS2,qS1,qS2,I) * dlkmC(l,dpS+dpI,M1-M2,0,psi,0) for l in [0,2]]) + \
                        sum([R_g(l,jK1,jK2,L1,L2,K1,K2) * w3jmatlabC(L1,l,L2,M1,M2-M1,-M2) * B0 * \
                             Ax_gC(l,pI1,pI2,qI1,qI2,pS1,pS2,qS1,qS2,I) * dlkmC(l,dpS+dpI,M1-M2,0,psi,0) for l in [0,2]]) \
                        )
#                mat[j,i] = mat[i,j]   #forcing symmetry in the matrix
    
    mat = mat.tocsr()
    return mat
