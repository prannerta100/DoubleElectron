from __future__ import division
from scipy.sparse import lil_matrix
from common_defs import *
def DiffTensorNoPotential(ind_arr, R, g):
    g0 = sum(g)/3
    cfac = g0 * betae / hbar
    #print(cfac)
    ndim = len(ind_arr)
    if len(ind_arr[0]) != 9: #13: #LKM jM jK p,q 2 electrons, 2 nuclei; 2x4+5=13
        raise "Problem with indices!"
    mat = lil_matrix((ndim,ndim), dtype=float)
    Rx, Ry, Rz = R[0]/cfac, R[1]/cfac, R[2]/cfac
    #print(Rx,Ry,Rz)
    for i in range(ndim):
        L = ind_arr[i][0] #first element
        K = ind_arr[i][2] #third element
        mat[i,i] = 0.5 * (Rx+Ry) * (L*(L+1)-K*K) + Rz * K*K
        if Rx != Ry:
            x = list(ind_arr[i])
            x[2] += 2
            if x in ind_arr:
                j = ind_arr.index(x)
                mat[i,j] = (Rx - Ry) * NormFacMinus(L,K+2)/(4*NormFacK(K,K+2))
                mat[j,i] = mat[i,j]
            #x = list(x0)
            #x[2] -= 2
            #if x in ind_arr:
            #    j = ind_arr.index(x)
            #    mat[i,j] = (Rx - Ry) * NormFacM(L,K-2)/(4*NormFacK(K,K-2))
    mat = mat.tocsr() #convert to csr format
    return mat
