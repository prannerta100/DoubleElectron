#bit of philosophy here: the goal is to just have a master template for how the starting vectors, vectors after pulse propagator, etc. will look
#the starting vector has the following hierarchy: L jK K M pIa qIa <pSa,qSa> pIb qIb <pSb,qSb>
#let us do just 1 3 pulse coherence pathway: +1=(+1,0;0,+1) ---> +2=(+1,+1) ----> -1=(0,-1;-1,0)
#we will do the following: stv is in +1=(+1,0;0,+1), so create a vector, define a matvec function that reshapes and multiplies, 
#the matvec function has 3 parts to it: electron 1 Hamiltonian, electron 2 Hamiltonian, dipolar term

#we will do the following: stv is in +1=(+1,0;0,+1), so create a vector, define a matvec function that reshapes and multiplies
#ideally I need 2 reshape functions, one for electron 1 and one for electron 2
#a nice benchmark would be to test this reshape and multiply functionality

#we then propagate to +1,+1, the pulse propagator for a pi/2 pulse should be the same as Lee 1994, then do matvec again, this time use pS=+1, pS=+1 matrices for both electrons
#then we need to convert +1,+1 to -1=(0,-1;-1,0), which is another pulse propagator

#there is some t_dip, etc. business to look into as well; so to be continued

from common_defs import *
def indgen(Lemax, Lomax, Kmax, Mmax, ipnmx, g, A, I2, angd, angm, angdip, psi):
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

    #dipolar tilt (all angles in degrees) !!!NEW!!!
    aldip = angdip[0]#-5
    bedip = angdip[1]#10
    gadip = angdip[2]#-15
    itdip = abs(aldip) > rndoff or abs(bedip) > rndoff or abs(gadip) > rndoff

    ipsi0 = (abs(psi)>rndoff) #director tilt flag; multiple ort or given director tilt
    kptmx = 0 #I set it, but I don't see it defined to be 0 anywhere. They define it only when c20, etc. are on. It is argmax_K c_{LK}.
    delta_K = 1 + (itd == 0 and itm == 0 and itdip == 0)   #!!!NEW!!!
    delta_L = 1 + (axialA and axialg and delta_K == 2 and kptmx == 0)
    #jMmin = 1 #as Zeeman Nuclear is 0
    jKmin = 1 - 2 * (abs(alm) > rndoff or abs(gam) > rndoff or abs(gad) > rndoff or abs(gadip) > rndoff) #as per rules from MIIMFreed, I read this in the code

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
#so everything before this line must be common to both electrons
#we will get into the loop, and then just add a subloop for the 2 other electrons

    ind_offdiag = []
    ind_diag = []
    ind_Plus1 = []
    ind_Plus2 = []
    ind_Minus1 = []

    ndimPlus1 = 0
    ndimPlus2 = 0
    ndimMinus1 = 0
    Lstarts_Plus1 = []
    Lstarts_Plus2 = []
    Lstarts_Minus1 = []
    
    ndimo = 0
    ndimd = 0
    Lstarts_offdiag = []
	Lstarts_diag = []

    Llist = [] #common, used for truncating the for loops over column indices
    #this loop comes from fbasis.f in nlspmc, just removed the jM parts to make it work for the jk basis
    #Question! Will this work for the dipolar superoperator? There we have the dipolar angle, i.e., angle between the dipolar and the lab frame  dipolar---> diffusion --> Lab; D2_00(dip--> Lab) = \Sum_m D2_0m(dip --> diff) D2_m0(diff--> L)
    for L in range(0,Lemax+1,delta_L):
        if par(L) == -1 and L > Lomax:
            continue
        Llist.append(L)
        Lstarts_offdiag.append(ndimo)
        Lstarts_diag.append(ndimd)
        Lstarts_Plus1.append(ndimPlus1)
        Lstarts_Plus2.append(ndimPlus2)
        Lstarts_Minus1.append(ndimMinus1)

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
                    for pIa in range(pImin, pImax+1):
                        #if M==0 and pI==0 and jM != par(L):
                        #    continue
                        if ipsi0==0 and pIa != M:
                            continue
                        qImaxa=I2 - abs(pIa)
                        qImina=-qImaxa
                        for qIa in range(qImina, qImaxa+1, 2):
                            ndimo += 1
                            #ind_diag, ind_offdiag used for merely creating the matrices
                            ind_offdiag.append([L,M,K,jK,pIa,qIa,1,0]) #jM removed
                            ndimd += 2
                            ind_diag.append([L,M,K,jK,pIa,qIa,0,1]) #jM removed
                            ind_diag.append([L,M,K,jK,pIa,qIa,0,-1]) #jM removed
                            for pIb in range(pImin, pImax+1):
                                #if M==0 and pI==0 and jM != par(L):
                                #    continue
                                if ipsi0==0 and pIa != M:
                                    continue
                                qImaxa=I2 - abs(pIa)
                                qImina=-qImaxa
                                for qIa in range(qImina, qImaxa+1, 2):
                                    #Plus1=(0,1),(1,0)
                                    ndimPlus1 += 4
                                    ind_Plus1.append([L,M,K,jK,pIa,qIa,1,0,pIb,qIb,0,1]) #jM removed
                                    ind_Plus1.append([L,M,K,jK,pIa,qIa,1,0,pIb,qIb,0,-1]) #jM removed
                                    ind_Plus1.append([L,M,K,jK,pIa,qIa,0,1,pIb,qIb,1,0]) #jM removed
                                    ind_Plus1.append([L,M,K,jK,pIa,qIa,0,-1,pIb,qIb,1,0]) #jM removed
                                    #Plus2=(1,1) 
                                    ndimPlus2 += 1
                                    ind_Plus2.append([L,M,K,jK,pIa,qIa,1,0,pIb,qIb,1,0]) #jM removed
                                    #Minus1=(0,-1),(-1,0)
                                    ndimMinus1 += 4
                                    ind_Minus1.append([L,M,K,jK,pIa,qIa,-1,0,pIb,qIb,0,1]) #jM removed
                                    ind_Minus1.append([L,M,K,jK,pIa,qIa,-1,0,pIb,qIb,0,-1]) #jM removed
                                    ind_Minus1.append([L,M,K,jK,pIa,qIa,0,1,pIb,qIb,-1,0]) #jM removed
                                    ind_Minus1.append([L,M,K,jK,pIa,qIa,0,-1,pIb,qIb,-1,0]) #jM removed

    #so understand this carefully: earlier L M K and L -M K used to show up next to each other, now that won't happen, 
    #now qS=-+1 are neighbors, changed the bandedness of the matrix a bit, anyway L,M,K are the hard hitters 
    #i.e., L1,M1,K1 can't be too far from L2,M2,K2	
    Lstarts_offdiag.append(ndimo) #caution: Lstarts lists have a length 1 + len(Llist), this is to take care of edge cases
    Lstarts_diag.append(ndimd)

    return ndimo, ndimd, ndimPlus1, ndimPlus2, ndimMinus1, \
           ind_offdiag, ind_diag, ind_Plus1, ind_Plus2, ind_Minus1, \
           Lstarts_offdiag, Lstarts_diag, Lstarts_Plus1, Lstarts_Plus2, Lstarts_Minus1,\
           Llist #(just a list of L's, doesn't matter whether it is single-electron or double-electron coherence) 
