# -*- coding: utf-8 -*-
"""
Created on Tue Sep 08 20:16:27 2015

@author: mike

# Cutoff for the solar atmosphere using the VALIIIc Model

##Background - the problem



	

Compute the atmospheric cut-off frequenccy (From Roberts 2004)

$$ \omega_{c}=\frac{\gamma g}{4\pi c_{s}}\sqrt{1+2\frac{d}{dz}\frac{P}{\rho g}}$$

The atmospheric scale parameter 

Read VAL IIc data from csv file

column 1 height [m]
Column 2 Temperature [K]
Column 3 Density [kg/m^3]
Column 4 Pressure [Pa]






"""

import numpy
import matplotlib.pyplot as plt
import math
from math import exp, sqrt, pi
#from math import exp, sqrt, pi



class vars:  
    mu=0.6e0
    R=8.31e3
    fgamma=1.66666667e0
    ggg=-274.0e0
    mu=4*pi/1.0e7
    

def cs(var, p, rho):
    vcs=sqrt(vars.fgamma*p/rho)
    return vcs

def lambda0(vars, P,rho):
    vlam0=P/(rho*vars.ggg)
    return vlam0

def lagrange_interp(xval,f,x,i): 
    t1=(xval-x[i])*(xval-x[i+1])/((x[i-1]-x[i])*(x[i-1]-x[i+1]))
    t2=(xval-x[i-1])*(xval-x[i+1])/((x[i]-x[i-1])*(x[i]-x[i+1]))
    t3=(xval-x[i-1])*(xval-x[i])/((x[i+1]-x[i-1])*(x[i+1]-x[i]))
    y=t1*f[i-1]+t2*f[i]+t3*f[i+1]
    return y

def diff5p(y,i,h):
    diff=(y[i-2]-8*y[i-1]+8*y[i+1]-y[i+2])/(12*h)
    return diff

def diff3p(y,i,h):
    diff=(y[i+1]-y[i-1])/(2*h)    
    return diff


#%matplotlib inline 
data=numpy.loadtxt(fname='..\\data\\atmos.csv', delimiter=',')
#print data
print data.shape
#print data[1:10,2]
height=data[0:2048,0]
ltemp=numpy.log10(data[0:2048,1]) #notation to avoid http://stackoverflow.com/questions/21687581/typeerror-only-length-1-arrays-can-be-converted-to-python-scalars-while-trying

print ltemp
plt.plot(height/1.0e6,ltemp)




ptemp=numpy.float64(data[1:2048,3])
rhotemp=numpy.float64(data[1:2048,2])

asize=data.size
print asize
acs = numpy.zeros(asize/4)
for i in range(0,(asize/4)-1):
    acs[i]=cs(vars,ptemp[i],rhotemp[i]);
    
print acs.size
print height.size
plt.plot(height/1.0e6,acs/1000)

#compute lambda0
alam0 = numpy.zeros(asize/4)
alami0 = numpy.zeros(asize/4)

for i in range(0,(asize/4)-1):
    alam0[i]=lambda0(vars,ptemp[i],rhotemp[i])

dh=height[0]-height[1]    
for i in range(1,(asize/4)-2):
    xval=height[0]-i*dh
    alami0[i]=lagrange_interp(xval,alam0,height,i)    
    
    
#print alam0  
#compute cutoff
atc0 = numpy.zeros(asize/4)
lamdash0 = numpy.zeros(asize/4)
for i in range(2,(asize/4)-4):
    h=height[i]-height[i+1]
    #lamdash0[i]=diff3p(alam0,i,h)
    #lamdash0[i]=ptemp[i]/rhotemp[i]
    lamdash0[i]=diff5p(alami0,i,h)
    #print h,alam0[i],lamdash0[i]
    print h,ptemp[i],rhotemp[i],lamdash0[i]
    atc0[i]=1.0/((vars.fgamma*vars.ggg/(4*pi*acs[i]))*math.sqrt(1+2*lamdash0[i]))
    #atc0[i]=1.0/((vars.fgamma*vars.ggg/(4*pi*acs[i])))
    #atc0[i]=math.sqrt(1+2*lamdash0[i])
    
#print atc0
plt.plot(height/1.0e6,lamdash0)