import numpy as np
import math
from LSymBasisModule import w3jmatlabC, dlkmC, Ax_dipC
from scipy.sparse import lil_matrix
from common_defs import *

def SecularDipolarTensor(ind_arr, D, angdip, Llist, Lstarts):
    #aldip,bedip,gadip a new set of Euler angles
    aldip = angdip[0]
    bedip = angdip[1]
    gadip = angdip[2]

    #massaging D (need to convert to Gauss?)
    D *= 1

    #  z   z       +  -     -  +                                     2 
    # S   S    - (S  S   + S  S  ) / 4     Hamiltonian, times D and D  (omega_L-->dip)
    #  1   2       1  2     1  2                                     00
    
    #get the spin part and the Wig 3j symbols; selection rules: I hope those from the previous diff (added potential?) and spin Hamiltonian operators extend here 
    #but definitely reconsider the selection rules before submitting the paper 
    #analogous to the Ax_a(2,0) term with nuclear spin I replaced by S2  
    #Dipolar term is D * Ax_a(2,0) * D2_00(omega_L --> dip)
    #       =D *    Ax_a(2,0)      * D2_0m(omega)            *                  D2_m0(angdip)
    #Mat el: D * <1|Ax_a(2,0)^x|2> * (L1 2  L2) * (M1==M2)   * (L1 2      L2) * D2_(K1-K2)0(angdip) * (-1)^(M1+K1) * N_L(L1,L2)
    #                                (M1 0 -M2)                (K1 K2-K1 -K2)
    #now the K symmetrization
    #Mat el: {Mat_el(K1,K2) + Mat_el(K1,-K2) jK2 (-1)^(L2+K2) + Mat_el(-K1,K2) jK1 (-1)^(L1+K1) + Mat_el(-K1,-K2) jK1jK2 (-1)^(L1+L2+K1+K2)}
    #               * N_K(K1,K2) * sqrt(jK1*jK2)
    #= D * <1|Ax_a(2,0)^x|2> * (L1 2  L2) * (M1==M2)   * (-1)^M1+K1 * N_L(L1,L2) * {

    #             (L1 2      L2) * D2_(K1-K2)0(angdip) + 3 more terms } no more simplification [might make execution time 2x/4x]
    #             (K1 K2-K1 -K2)
    ndim = len(ind_arr)
    if len(ind_arr[0]) != 12:
        raise NameError("Problem with indices!") #just a check
    mat = lil_matrix((ndim,ndim),dtype=complex) #complex because of something, not sure what 

    for i in range(ndim): #ndim is ndimd
        L1,M1,K1,jK1,pIa1,qIa1,pSa1,qSa1,pIb1,qIb1,pSb1,qSb1 = ind_arr[i]
        #print(L1,M1,K1,jK1,pI1,qI1,pS1,qS1)
        #if pS1 != 0 or qS1 != 1:
        #    raise ValueError("wrong routine used, check index list again")
        #left limit
        if L1-2 in Llist:
            left_j = Lstarts[Llist.index(L1-2)]
        elif L1-1 in Llist:
            left_j = Lstarts[Llist.index(L1-1)]
        else:
            left_j = Lstarts[Llist.index(L1)]
        #right limit
        if L1+2 in Llist:
            right_j = Lstarts[Llist.index(L1+2)+1]#note the starting right after L1+2
        elif L1+1 in Llist:
            right_j = Lstarts[Llist.index(L1+1)+1]#Lstarts arrays have a length 1+len(Llist), last element must be ndim
        else:
            right_j = ndim
        #print('l,r,n',left_j,right_j,ndim) 
        for j in range(left_j,right_j):  #range(ndim) #define band, this will really speed up things! Think about the nbnd in CW nlsl matrll.f
            L2,M2,K2,jK2,pIa2,qIa2,pSa2,qSa2,pIb2,qIb2,pSb2,qSb2 = ind_arr[j]
            if M1==M2 and pIa1==pIa2 and pIb1==pIb2 and qIa1==qIa2 and qIb1==qIb2 and abs(K2-K1) <=2:
                #Mat el: D * <1|Ax_a(2,0)^x|2> * (L1 2  L2) * (M1==M2)   * (L1 2      L2) * D2_(K1-K2)0(angdip) * (-1)^(M1+K1) * N_L(L1,L2)
                #                                (M1 0 -M2)                (K1 K2-K1 -K2)
                mat[i,j] = D * NormFacL(L1,L2) * 0.5 * NormFacK(K1,K2) * par(M1) * w3jmatlabC(L1,2,L2,M1,0,-M2) *\
                           Ax_dipC(pSa1, pSa2, qSa1, qSa2,  pSb1, pSb2, qSb1, qSb2) * math.sqrt(jK1*jK2) * (\
                                                  dlkmC(2,K1-K2,0,aldip,bedip,gadip)*par(K1)*w3jmatlabC(L1,2,L2,K1,K2-K1,-K2) + \
                                             jK1* dlkmC(2,-K1-K2,0,aldip,bedip,gadip)*par(L1)*w3jmatlabC(L1,2,L2,-K1,K2+K1,-K2) + \
                                             jK2* dlkmC(2,K1+K2,0,aldip,bedip,gadip)*par(L2+K2+K1)*w3jmatlabC(L1,2,L2,K1,-K2-K1,K2) + \
                                         jK1*jK2* dlkmC(2,-K1+K2,0,aldip,bedip,gadip)*par(L2+L1+K2)*w3jmatlabC(L1,2,L2,-K1,-K2+K1,K2) )
                                
                #using Eva's 1982 paper
#                mat[j,i] = mat[i,j] #forced symmetrization, shouldn't have to use it

    mat = mat.tocsr()
    return mat

