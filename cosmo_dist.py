#THIS CODE IS INDENTED TO CALCULATE COSMOLOGICAL DISTANCES 
#UPDATE: 6 JULY 2015

from math import *
import numpy as np
from scipy import integrate
import sys

#INPUTS: H, 0_M, O_L, z (Hubble constant, omega of mass, omega of lambda and the redshift
try:
    if sys.argv[1] == '-h':
        print '''Cosmological calculator (based off Hogg 2000)
              input values =  Ho, redshift, Omega_m, Omega_lambda
              ouput values = age at z, distance in Mpc, kpc/arcsec, apparent to abs mag conversion
              
              Options:
              -h for this message 
              -v for verbose response '''
        sys.exit()
    if sys.argv[1] == '-v':
        verbose = 1
        length = len(sys.argv)-1
    else:
        verbose = 0
        length = len(sys.argv)
    if (length != 5):
        print 'need some values or too many values'
        sys.exit()
    elif (length == 5):
        H = float(sys.argv[1+verbose])      #Hubble Constant (km/s/Mpc)
        z =  float(sys.argv[2+verbose])     #redshift
        O_M = float(sys.argv[3+verbose])    #omega of mass
        O_L = float(sys.argv[4+verbose])   #omega of lambda 


    O_K = 1.0 - O_M - O_L #curvature of space
    c = 299792.458       #speed of Light (Km/s)
    D_H = c/H        #Hubble Distance 

    #COMOVING DISTANCE (Mpc) (Credit: Hogg 2000 Eq15)
    integ = lambda z: (np.sqrt(O_M*(1.0+z)**3+O_K*(1.0+z)**2+O_L))**-1
    D_C, err = integrate.quad(integ,0,z)
    D_C = D_H*D_C
    
    #other way to do the integrate
    N = 200
    dz = z/N
    z_k = 0.0
    z_k1 = z_k+dz
    integrad = 0.0
    while (z_k1 <= z+dz):
        integrad = integrad + (np.sqrt(O_M*(1.0+z_k)**3+O_K*(1.0+z_k)**2+O_L))**-1 + (np.sqrt(O_M*(1.0+z_k1)**3+O_K*(1.0+z_k1)**2+O_L))**-1
        z_k = z_k+dz
        z_k1 = z_k1+dz
        
    test = (z/(2.0*N))*integrad
    print ""
    print "the difference between the python integral function and the trapezoidal rule (for calculating the comoving distance) is:", D_C-test*D_H, "Mpc"

    #TRANSVERSE COMOVING DISTANCE (Mpc) (Credit: Hogg 2000 Eq15)
    if (O_K > 0.0):
        D_M = (D_H/np.sqrt(O_K))*np.sinh(np.sqrt(O_K)*D_C/D_H)

    elif (O_K == 0.0):
        D_M = D_C
    
    elif (O_K < 0.0):
        D_M = (D_H/np.sqrt(np.absolute(O_K)))*np.sin(np.sqrt(np.absolute(O_K))*D_C/D_H)

    #ANGULAR DIAMETER DISTANCE (Credit: Hogg 2000 Eq18)
    D_A = D_M / (1.0+z)

    #LUMINOSITY DISTANCE (Credit: Hogg 2000 Eq15)
    D_L = (1.0+z)*D_M

    #DISTANCE MODULUS (pc) (Credit: Hogg 2000 Eq25)
    DM = 5.0*np.log10((D_L*(1.0E6))/10.0)


    #LIGHT TRAVEL DISTANCE (Credit: Hogg 2000 Eq30)
    integ = lambda z: ((1.0+z)*((np.sqrt(O_M*(1.0+z)**3+O_K*(1.0+z)**2+O_L))))**-1
    D_T, err = integrate.quad(integ,0,z)
    D_T = D_H*D_T

    #other way to do the integrate
    N = 200
    dz = z/N
    z_k = 0.0
    z_k1 = z_k+dz
    integrad = 0.0
    while (z_k1 <= z+dz):
        integrad = integrad + ((1.0+z_k)*((np.sqrt(O_M*(1.0+z_k)**3+O_K*(1.0+z_k)**2+O_L))))**-1 + ((1.0+z_k1)*((np.sqrt(O_M*(1.0+z_k1)**3+O_K*(1.0+z_k1)**2+O_L))))**-1
        z_k = z_k+dz
        z_k1 = z_k1+dz
        
    test = (z/(2.0*N))*integrad
    print ""
    print "the difference between the python integral function and the trapezoidal rule (for calculating the light travel distance) is:", D_T-test*D_H, "Mpc"
    
    #COMOVING VOLUME (Credit: Hogg 2000 Eq29)
    if (O_K > 0.0):
        V_C = ( (4.0*np.pi*D_H**3)/(2.0*O_K) )*( (D_M/D_H)*np.sqrt(1.0+O_K*(D_M**2/D_H**2))-(1.0/np.sqrt(np.absolute(O_K)))*np.sinh(np.absolute(O_K)*(D_M/D_H)))

    elif (O_K == 0.0):
        V_C = ((4.0*np.pi)/3.0)*D_M**3

    elif (O_K < 0.0):
        V_C = ( (4.0*np.pi*D_H**3)/(2.0*O_K) )*( (D_M/D_H)*np.sqrt(1.0+O_K*(D_M**2/D_H**2))-(1.0/np.sqrt(np.absolute(O_K)))*np.sin(np.absolute(O_K)*(D_M/D_H)))

    #Printing all the outputs
    print ""
    print "the transverse comoving distance is",D_M,"Mpc"
    print "the angular diameter distance is", D_A,"Mpc"
    print "the luminosity distance is", D_L,"Mpc"
    print "the distance modulus is", DM
    print "the comoving volume is", V_C, "Mpc^3 or", V_C*1e-9,"Gpc^3"
    print ""
    print "D_M/D_H:",D_M/D_H
    print "D_A/D_H", D_A/D_H
    print "D_L/D_H", D_L/D_H
    print "D_T/D_H", D_T/D_H
except IndexError:
    print 'need some values or too many values'
