#!/usr/bin/env python
# Plots data in binary dumps from sac smaug.
# Usage: python sac-rea.py ../../../hermes-temp/bin/nl-shwave.0024.bin

from numpy import *
import sys

from matplotlib import use
#use('Agg') 
from matplotlib import pyplot
from pylab import *
import struct

#
# Read binary file
#
try:
  file = open(sys.argv[1],'rb')
except:
  print 'Usage: ./sac-read.py <binary_dump>'
  raise SystemExit

file.seek(0,2)
eof = file.tell()
file.seek(0,0)
#name=file.readline()
#buffer1=file.read(4)
#file.read.string()
#hr1=file.read(4)
name = file.read(79)
#hr2=file.read(4)

#hr1=file.read(4)
nit = fromfile(file,dtype=int32,count=1)
#buffer2=file.read(4)
#l = struct.unpack( '>i',buffer2 )
#buffer3=file.read(8)
#t = struct.unpack( '>d',buffer3 )
#n=fromstring(buffer1, dtype=int64,count=1)
#t=fromstring(buffer2, dtype=double,count=1)
t = fromfile(file,dtype=float64,count=1)
ndim=fromfile(file,dtype=int32,count=1)
neqpar=fromfile(file,dtype=int32,count=1)
nw=fromfile(file,dtype=int32,count=1)
#hr1=file.read(4)
#hr2=file.read(4)
#nit = fromfile(file,dtype=int32,count=1)

#nx = struct.unpack( '>i',file.read(4) )
#ny = struct.unpack( '>i',file.read(4) )
ndata = fromfile(file,dtype=int32,count=ndim)[:ndim]
#nx= ndata[0]
#ny= ndata[1]
#nz= ndata[2]
varbuf = fromfile(file,dtype=float,count=7)[:7]

#if ndim=2
varnames = file.read(79)

#if ndim=3
if ndim==3:
    alldat=fromfile(file,dtype=float,count=(nw+ndim)*ndata[0]*ndata[1]*ndata[2])[:(nw+ndim)*ndata[0]*ndata[1]*ndata[2]]
    #if size(alldat)<(nw+ndim)*ndata[0]*ndata[1]*ndata[2]:
    #    alldat=resize(alldat,(nw+ndim)*ndata[0]*ndata[1]*ndata[2])
    alldat=np.reshape(alldat,(nw+ndim,ndata[0],ndata[1],ndata[2],),'C')
elif ndim==2:
    alldat=fromfile(file,dtype=float,count=(nw+ndim)*ndata[0]*ndata[1])[:(nw+ndim)*ndata[0]*ndata[1]]
    if size(alldat)<(nw+ndim)*ndata[0]*ndata[1]:
        alldat=resize(alldat,(nw+ndim)*ndata[0]*ndata[1])
    alldat=np.reshape(alldat,(nw+ndim,ndata[0],ndata[1]))







#if file.tell() != eof: print 'Error: Too few bytes read.'

file.close()


#display config information
print name
#print shape
#print nx,ny,nz
#print ndata
#print t,dt
#print gamma1,cs
#print coordsys
#print count


#print y[0],y[1],y[191]
#print rho.shape
#print rho[0][0],rho[191][191]

