# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 16:58:21 2019

@author: ppg22
"""
from __future__ import division
import numpy as np
from sympy.physics.quantum.spin import Rotation
from sympy import N
def dlkm_orig(l,k,m,alpha,beta,gamma):
    #alpha, beta, gamma in degrees
    if l not in [0,2]: 
        return 0.0
    if l==0:
        return 1.0
    if abs(k)>l or abs(m)>l:
        return 0
    if l==2:
        alpha *= np.pi/180.0
        beta *= np.pi/180.0
        gamma *= np.pi/180.0
        #print(l,k,m,"func_print_from_dlkm")
        ans = complex(N(Rotation.D(l,k,m,alpha,beta,gamma).doit()))
        return ans

def dlkm(l,k,m,alpha,beta,gamma):
    #alpha, beta, gamma in degrees
    if l not in [0,2]: 
        return 0.0
    if abs(k)>l or abs(m)>l:
        return 0
    if l==0:
        return 1.0
    if l==2:
        alpha *= np.pi/180.0
        beta *= np.pi/180.0
        gamma *= np.pi/180.0
        cb = np.cos(beta)
        sb = np.sin(beta)
        #print(l,k,m,"func_print_from_dlkm")
        if (k==2 and m==2) or (k==-2 and m==-2):
            d = (1+cb) * (1+cb) / 4.0
        elif (k==2 and m==1) or (k==-1 and m==-2):
            d = -sb * (1+cb) / 2.0
        elif (k==1 and m==2) or (k==-2 and m==-1):
            d = sb * (1+cb) / 2.0
        elif (k==-2 and m==0) or (k==0 and m==-2) or (k==2 and m==0) or (k==0 and m==2):
            d = np.sqrt(3.0/8.0) * sb * sb
        elif (k==2 and m==-1) or (k==1 and m==-2):
            d = sb * (cb - 1) / 2.0
        elif (k==-2 and m==1) or (k==-1 and m==2):
            d = -sb * (cb - 1) / 2.0
        elif (k==2 and m==-2) or (k==-2 and m==2):
            d = (1-cb) * (1-cb) / 4.0
        elif (k==1 and m==1) or (k==-1 and m==-1):
            d = (2*cb-1) * (1+cb) / 2.0
        elif (k==1 and m==-1) or (k==-1 and m==1):
            d = (2*cb+1) * (1-cb) / 2.0
        elif (k==1 and m==0) or (k==0 and m==-1):
            d = -sb * cb * np.sqrt(1.5)
        elif (k==0 and m==1) or (k==-1 and m==0):
            d = sb * cb * np.sqrt(1.5)
        elif k==0 and m==0:
            d = (3*cb*cb - 1) / 2.0
        ans = d * complex(np.cos(k * alpha + m * gamma), -np.sin(k * alpha + m * gamma)) 
        return ans

