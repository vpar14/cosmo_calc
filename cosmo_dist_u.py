#THIS CODE IS INDENTED TO CALCULATE COSMOLOGICAL DISTANCES 
#UPDATE: 6 JULY 2015
from math import *
import numpy as np
from scipy import integrate
import sys

class Cos_Para(object):
    def __init__(self, H=70.2,z=0,O_M=0.725,O_L=0.275):
        self.H = H
        self.z = z
        self.O_M=O_M
        self.O_L=O_L
    
    def O_K(self):
        return 1.0-self.O_M-self.O_L

    def D_H(self):
        return 299792.458/self.H

    def integ(self):
        return lambda z:(np.sqrt(self.O_M*(1.0+self.z)**3+self.O_K*(1.0+self.z)**2+self.O_L))**-1

    def D_C(self):
        #COMOVING DISTANCE (Mpc) (Credit: Hogg 2000 Eq15)
        return self.D_H*integrate.quad(self.integ,0,self.z)

    def D_M(self):
        if self.O_K > 0.0:
            return (self.D_H/np.sqrt(self.O_K))*np.sinh(np.sqrt(self.O_K)*self.D_C/self.D_H)
        elif self.O_K == 0.0:
            return self.D_C
        elif self.O_K < 0.0:
            return (self.D_H/np.sqrt(np.absolute(self.O_K)))*np.sin(np.sqrt(np.absolute(self.O_K))*self.D_C/self.D_H)

    def D_A(self):
        return self.D_M / (1.0+self.z)

    def D_L(self):
        return (1.0+self.z)*self.D_M

    def distance_modulus(self):
        return 5.0*np.log10((self.D_L*(1.0E6))/10.0) 

    def D_T(self):
        return self.D_H * integrate.quad(self.integ/(1+self.z),0,self.z)

    def V_C(self):
        if self.O_K > 0.0:
            return ( (4.0*np.pi*self.D_H**3)/(2.0*self.O_K) )*( (self.D_M/self.D_H)*np.sqrt(1.0+self.O_K*(self.D_M**2/self.D_H**2))-(1.0/np.sqrt(np.absolute(self.O_K)))*np.sinh(np.absolute(self.O_K)*(self.D_M/self.D_H)))
        elif (self.O_K == 0.0):
            return ((4.0*np.pi)/3.0)*self.D_M**3 
        elif (self.O_K < 0.0):
            return ( (4.0*np.pi*self.D_H**3)/(2.0*self.O_K) )*( (self.D_M/self.D_H)*np.sqrt(1.0+self.O_K*(self.D_M**2/self.D_H**2))-(1.0/np.sqrt(np.absolute(self.O_K)))*np.sin(np.absolute(self.O_K)*(self.D_M/self.D_H))) 

