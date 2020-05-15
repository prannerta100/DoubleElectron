import numpy as np
import math
from scipy.sparse import coo_matrix
from scipy.sparse import kron as spkron
from scipy.sparse import eye as speye
def PulsePropPlus1ToPlus2(ndimPlus1, ndimPlus2, theta=math.pi/2, phi=0.0):
    prefac =  0.5j * math.cos(theta/2)**2 * math.sin(theta) * (math.cos(phi) - 1.0j*math.sin(phi))
    #prefac: i/2   * cos^2(theta/2)       *      sin(theta) *  exp(-i \phi) -> pg 31 lab nb 04/20
    if ndimPlus1 != ndimPlus2 * 4:
        raise ValueError('ndimPlus1 != 4 x ndimPlus2') #caution!
    I = np.kron(np.arange(ndimPlus2),np.ones(4)).astype(int)
    J = np.arange(ndimPlus1).astype(int)
    E = prefac * np.kron(np.ones(ndimPlus2),np.array([1.0,-1.0,1.0,-1.0])).astype(complex)
    return coo_matrix((E,(I,J)), shape=(ndimPlus2,ndimPlus1))
    
def PulsePropPlus2ToMinus1(ndimPlus2, ndimMinus1, theta=math.pi/2, phi=0.0):
    prefac = 0.5j * math.sin(theta/2)**2 * math.sin(theta) * (math.cos(3*phi) + 1.0j*math.sin(3*phi))
    #prefac: i/2   * sin^2(theta/2)       *      sin(theta) *  exp(3i \phi) -> pg 31 lab nb 04/20
    if ndimMinus1 != ndimPlus2 * 4:
        raise ValueError('ndimMinus1 != 4 x ndimPlus2') #caution!
    I = np.arange(ndimMinus1).astype(int)
    J = np.kron(np.arange(ndimPlus2),np.ones(4)).astype(int)
    E = prefac * np.kron(np.ones(ndimPlus2),np.array([1.0,-1.0,1.0,-1.0])).astype(complex)
    return coo_matrix((E,(I,J)), shape=(ndimMinus1,ndimPlus2))

def PulsePropPlus2ToMinus2(ndimPlus2, ndimMinus2, theta=math.pi/2, phi=0.0):
    prefac = math.sin(theta/2)**4 * (math.cos(4*phi) + 1.0j*math.sin(4*phi))
    #prefac:      sin^2(theta/2)^4        exp(4i \phi) -> pg 35 lab nb 05/20
    if ndimPlus2 != ndimMinus2:
        raise ValueError('ndimPlus2 != ndimMinus2') #caution!
    I = np.arange(ndimPlus2).astype(int)
    J = np.arange(ndimMinus2).astype(int)
    E = prefac * np.ones(ndimPlus2).astype(complex)
    return coo_matrix((E,(I,J)), shape=(ndimMinus2,ndimPlus2))

def PulsePropMinus2ToMinus1(ndimMinus2, ndimMinus1, theta=math.pi/2, phi=0.0):
    prefac = 0.5j * math.cos(theta/2)**2 * math.sin(theta) * (math.cos(phi) - 1.0j*math.sin(phi))
    #prefac: i/2  *      cos^2(theta/2)  *      sin(theta) *  exp(-i \phi) -> pg 35 lab nb 05/20
    if 4*ndimMinus2 != ndimMinus1:
        raise ValueError('4 x ndimMinus2 != ndimMinus1') #caution!
    I = np.arange(ndimMinus1).astype(int)
    J = np.kron(np.arange(ndimMinus2),np.ones(4)).astype(int)
    E = prefac * np.kron(np.ones(ndimPlus2),np.array([-1.0,1.0,-1.0,1.0])).astype(complex)
    return coo_matrix((E,(I,J)), shape=(ndimMinus1,ndimMinus2))

def PulsePropPlus2ToPlus1(ndimPlus2, ndimPlus1, theta=math.pi/2, phi=0.0):
    prefac = 0.5j * math.cos(theta/2)**2 * math.sin(theta) * (math.cos(phi) + 1.0j*math.sin(phi))
    #prefac: i/2   * cos^2(theta/2)       *      sin(theta) *  exp(i \phi) -> pg 35 lab nb 05/20
    if ndimPlus1 != ndimPlus2 * 4:
        raise ValueError('ndimPlus1 != 4 x ndimPlus2') #caution!
    I = np.arange(ndimPlus1).astype(int)
    J = np.kron(np.arange(ndimPlus2),np.ones(4)).astype(int)
    E = prefac * np.kron(np.ones(ndimPlus2),np.array([1.0,-1.0,1.0,-1.0])).astype(complex)
    return coo_matrix((E,(I,J)), shape=(ndimPlus1,ndimPlus2))

def PulsePropPlus1ToMinus1(ndimPlus1, ndimMinus1, theta=math.pi/2, phi=0.0):
    prefac1 = 0.25 * math.sin(theta)**2 * (math.cos(2*phi) + 1.0j*math.sin(2*phi))
    prefac2 = math.sin(theta/2)**4 * (math.cos(2*phi) + 1.0j*math.sin(2*phi)) #not much common
    if ndimPlus1 != ndimMinus1:
        raise ValueError('ndimPlus1 != ndimMinus1') #caution!
    #to be continued, this part is tricky as there are 16 combinations
    mat4by4 = np.array([[prefac1,prefac2,prefac1,-prefac1],[prefac2,prefac1,-prefac1,prefac1],[prefac1,-prefac1,prefac1,prefac2],[-prefac1,prefac1,prefac2,prefac1]])
    return spkron(speye(int(ndimPlus1/4)),mat4by4,format='coo')
    
