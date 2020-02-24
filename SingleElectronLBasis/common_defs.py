import math
#rndoff
rndoff = 1e-15

#parity
def par(x):
    return 1 - 2 * (x % 2) 
#fundamental constants
hbar = 1.05443e-27
betae = 9.2731e-21

#normalization constants
def NormFacPlus(l,k):
    return math.sqrt((l-k-1)*(l-k)*(l+k+1)*(l+k+2))
def NormFacMinus(l,k):
    return math.sqrt((l+k-1)*(l+k)*(l-k+1)*(l-k+2))
def NormFacK(k1,k2):
    return 1/math.sqrt(((k1==0)+1)*((k2==0)+1))
def NormFacL(l1,l2):
    return math.sqrt((2*l1+1)*(2*l2+1))
def NormFacP(pI1,M1,pI2,M2):
    return 1/math.sqrt((1+(pI1==0)*(M1==0))*(1+(pI2==0)*(M2==0)))

