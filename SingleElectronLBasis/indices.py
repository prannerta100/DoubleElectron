from common_defs import *
def indgen(Lemax, Lomax, Kmax, Mmax, ipnmx, g, A, I2, angd, angm, psi):
    axialA = (A[0] == A[1]) #check if g axial
    axialg = (g[0] == g[1]) #check if A axial
    
    #diffusion tilt (all angles in degrees)
    ald = angd[0]#-5
    bed = angd[1]#10
    gad = angd[2]#-15
    itd = abs(ald) > rndoff or abs(bed) > rndoff or abs(gad) > rndoff

    #magnetic tilt
    alm = angm[0]#-20
    bem = angm[1]#25
    gam = angm[2]#-30
    itm = abs(alm) > rndoff or abs(bem) > rndoff or abs(gam) > rndoff

    ipsi0 = (abs(psi)>rndoff) #director tilt flag; multiple ort or given director tilt
    kptmx = 0 #I set it, but I don't see it defined to be 0 anywhere. They define it only when c20, etc. are on. It is argmax_K c_{LK}.
    delta_K = 1 + (itd == 0 and itm == 0)
    delta_L = 1 + (axialA and axialg and delta_K == 2 and kptmx == 0)
    #jMmin = 1 #as Zeeman Nuclear is 0
    jKmin = 1 - 2 * (abs(alm) > rndoff or abs(gam) > rndoff or abs(gad) > rndoff) #as per rules from MIIMFreed, I read this in the code

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
    Lstarts_offdiag = []
    Lstarts_diag = []
    Llist = []
    #this loop comes from fbasis.f in nlspmc
    for L in range(0,Lemax+1,delta_L):
        if par(L) == -1 and L > Lomax:
            continue
        Llist.append(L)
        Lstarts_offdiag.append(ndimo)
        Lstarts_diag.append(ndimd)
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

    return ndimo, ndimd, ind_offdiag, ind_diag, Lstarts_offdiag, Lstarts_diag, Llist
