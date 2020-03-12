import math
import numpy as np
from common_defs import *
from myModule import w3jmatlabC, dlkmC, Ax_aC

if abs(NormFacL(3,5) - math.sqrt(77)) > rndoff:
    raise ValueError("NormFacL unit test failed")
if abs(NormFacK(0,0) - 0.5) > rndoff or abs(NormFacK(0,2) - 1/math.sqrt(2)) > rndoff or abs(NormFacK(3,0) - 1/math.sqrt(2)) > rndoff or abs(NormFacK(4,5) - 1) > rndoff:
    raise ValueError("NormFacK unit test failed")

ald,bed,gad,alm,bem,gam, g, A = 5,10,15, 20,25,30, [complex(0,6),7,complex(36,71)], [complex(0,6),7,complex(36,71)]
d2diff = {m1:{m2:dlkmC(2,m1,m2,ald,bed,gad) for m2 in range(-2,3)} for m1 in range(-2,3)}
d2mag = {m1:{m2:dlkmC(2,m1,m2,alm,bem,gam) for m2 in range(-2,3)} for m1 in range(-2,3)}
d2diffmag = {m1:{m2:sum([d2diff[m1][mtmp] * d2mag[mtmp][m2] for mtmp in range(-2,3)]) for m2 in range(-2,3)} for m1 in range(-2,3)}
f2_aa = {0:math.sqrt(2/3)*(A[2]-0.5*(A[0]+A[1])), -1:0, 1:0, -2:0.5*(A[0]-A[1]), 2:0.5*(A[0]-A[1])}
f2_gg = {0:math.sqrt(2/3)*(g[2]-0.5*(g[0]+g[1])), -1:0, 1:0, -2:(g[0]-g[1])/2, 2:(g[0]-g[1])/2}
F_aD_Conj0 = np.conj(-sum(A)/math.sqrt(3))
F_gD_Conj0 = np.conj(-sum(g)/math.sqrt(3))
F_aD_Conj2 = {m:sum([d2diffmag[m][mtmp]* f2_aa[mtmp].conjugate() for mtmp in range(-2,3)]) for m in range(-2,3)}
F_gD_Conj2 = {m:sum([d2diff[m][mtmp] * f2_gg[mtmp].conjugate() for mtmp in range(-2,3)]) for m in range(-2,3)}

if abs(F_gD_Conj0 - -(43-77j)/math.sqrt(3)) > rndoff or abs(F_gD_Conj2[1]-(-5.311264597718748+11.996830758823908j)) > 1e-10:
    raise ValueError("F_gD_Conj unit test failed")
if abs(F_aD_Conj0 - -(43-77j)/math.sqrt(3)) > rndoff or abs(F_aD_Conj2[1]-(2.50766627945+37.0389445955j)) > 1e-10:
    raise ValueError("F_aD_Conj unit test failed")
#if abs(G_a(0,1,1,0) - -43/math.sqrt(3)) > rndoff or abs(G_a(0,-1,4,0) - 77/math.sqrt(3)) > rndoff:
#raise ValueError("G_a unit test failed")
#if abs(G_g(0,1,1,0) - -43/math.sqrt(3)) > rndoff or abs(G_g(0,-1,4,0) - 77/math.sqrt(3)) > rndoff:
#raise ValueError("G_g unit test failed")
#if abs(Ax_g(2,1,1,0,0)-0) > rndoff:
#raise ValueError("Ax_g unit test failed")
if abs(Ax_aC(2,1,0,0,1,1) - 3/8) > rndoff: #nuclear spin flip happening
    raise ValueError("Ax_a unit test failed")

#managed to reach here, means all tests passed!
print("Great, all unit tests passed")

