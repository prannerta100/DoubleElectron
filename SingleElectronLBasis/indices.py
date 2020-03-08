def indgen(Lemax, Lomax, Kmax, Mmax, ipnmx):
    ind_offdiag = []
    ind_diag = []
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
    pImaxval = ipnmx
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
#                for jM in range(jMmin,1+1,2):
                for M in range(-min(Mmax, L),min(Mmax, L)+1):
                    #calculate pImin and max
                    pImax=min(I2,pImaxval)
                    if M == 0:
                        pImin=0
                    else:
                        pImin=-pImax
                    for pI in range(pImin, pImax+1):
                        #if M==0 and pI==0 and jM != par(L):
                        #    continue
                        if ipsi0==0 and pI != M:
                            continue
                        qImax=I2 - abs(pI)
                        qImin=-qImax
                        for qI in range(qImin, qImax+1, 2):
                            ndimo += 1
                            ind_offdiag.append([L,M,K,jK,pI,qI,1,0]) #jM removed
                            ndimd += 2
                            ind_diag.append([L,M,K,jK,pI,qI,0,1]) #jM removed
                            ind_diag.append([L,M,K,jK,pI,qI,0,-1]) #jM removed
     #so understand this carefully: earlier L M K and L -M K used to show up next to each other, now that won't happen, 
     #now qS=-+1 are neighbors, changed the bandedness of the matrix a bit, anyway L,M,K are the hard hitters 
     #i.e., L1,M1,K1 can't be too far from L2,M2,K2	

    return ndimo, ndimd, ind_offdiag, ind_diag
