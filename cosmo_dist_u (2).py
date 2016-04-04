#THIS CODE IS INDENTED TO CALCULATE COSMOLOGICAL DISTANCES 
#UPDATE: 6 JULY 2015
from math import *
import numpy as np
from scipy import integrate
import sys

def integ(z,O_M,O_K,O_L):
    return (np.sqrt(O_M*(1.0+z)**3+O_K*(1.0+z)**2+O_L))**-1

def integ_2(z,O_M,O_K,O_L):
    return integ(O_M,O_K,O_L,z)/(1+z)

class Cos_Para(object):
    def __init__(self, H=70.2,z=0,O_M=0.275,O_L=0.725):
        self.H = H
        self.z = z
        self.O_M=O_M
        self.O_L=O_L
    
    @property
    def O_K(self):
        return 1.0-self.O_M-self.O_L

    @property
    def D_H(self):
        return 299792.458/self.H

    @property
    def D_C(self):
        return self.D_H*np.array(integrate.quad(integ,0,self.z,args=(self.O_M,self.O_K,self.O_L)))

    @property
    def D_M(self):
        if self.O_K > 0.0:
            return (self.D_H/np.sqrt(self.O_K))*np.sinh(np.sqrt(self.O_K)*self.D_C/self.D_H)
        elif self.O_K == 0.0:
            return self.D_C
        elif self.O_K < 0.0:
            return (self.D_H/np.sqrt(np.absolute(self.O_K)))*np.sin(np.sqrt(np.absolute(self.O_K))*self.D_C/self.D_H)

    @property
    def D_A(self):
        return self.D_M / (1.0+self.z)

    @property
    def D_L(self):
        return (1.0+self.z)*self.D_M

    @property
    def distance_modulus(self):
        return 5.0*np.log10((self.D_L*(1.0E6))/10.0) 

    @property
    def D_T(self):
        return self.D_H*np.array(integrate.quad(integ_2,0,self.z,args=(self.O_M,self.O_K,self.O_L)))

    @property
    def V_C(self):
        if self.O_K > 0.0:
            return ( (4.0*np.pi*self.D_H**3)/(2.0*self.O_K) )*( (self.D_M/self.D_H)*np.sqrt(1.0+self.O_K*(self.D_M**2/self.D_H**2))-(1.0/np.sqrt(np.absolute(self.O_K)))*np.sinh(np.absolute(self.O_K)*(self.D_M/self.D_H)))
        elif (self.O_K == 0.0):
            return ((4.0*np.pi)/3.0)*self.D_M**3 
        elif (self.O_K < 0.0):
            return ( (4.0*np.pi*self.D_H**3)/(2.0*self.O_K) )*( (self.D_M/self.D_H)*np.sqrt(1.0+self.O_K*(self.D_M**2/self.D_H**2))-(1.0/np.sqrt(np.absolute(self.O_K)))*np.sin(np.absolute(self.O_K)*(self.D_M/self.D_H))) 

