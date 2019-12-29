from sympy import N
from sympy.physics.wigner import *
from myModule import w3jmatlabC
print(N(wigner_3j(10,1,10,7,1,-8))- w3jmatlabC(10,1,10,7,1,-8))
print(N(wigner_3j(8,2,10,8,2,-10)) - w3jmatlabC(8,2,10,8,2,-10))
print(N(wigner_3j(15,12,21,0,1,-1))- w3jmatlabC(15,12,21,0,1,-1))
print(N(wigner_3j(2,2,1,1,-1,0)) - w3jmatlabC(2,2,1,1,-1,0))
print(N(wigner_3j(2,2,1,2,-2,0)) - w3jmatlabC(2,2,1,2,-2,0))
print(N(wigner_3j(4,2,3,2,-2,0)) - w3jmatlabC(4,2,3,2,-2,0))
print(N(wigner_3j(8,2,10,2,-1,-1)) - w3jmatlabC(8,2,10,2,-1,-1))
print(N(wigner_3j(1,1,2,0,1,-1)) - w3jmatlabC(1,1,2,0,1,-1))
