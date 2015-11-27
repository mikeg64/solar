from math import exp, sqrt, pi

class vars:
  
    mu=0.6e0
    R=8.31e3
    fgamma=1.66666667e0
    ggg=274.0e0
    mu=4*pi/1.0e7
    

def cs(var, p, rho):
    vcs=sqrt(vars.fgamma*p/rho)
    return vcs

ptemp=numpy.float64(data[1:2048,3])
rhotemp=numpy.float64(data[1:2048,2])
#ptemp=numpy.float64(301.2741)

asize=data.size
print asize
acs = numpy.zeros(asize/4)
for i in range(0,(asize/4)-3):
    #acs[i]=cs(vars,ptemp[i],rhotemp[i]);
    count=1
    ptot=ptemp[i]
    rhottot=rhotemp[i]
    for j in range(i+1,i+20):        
        if j<=(asize/4)-4:
            ptot=ptot+ptemp[j]
            rhottot=rhottot+rhotemp[j]
            count=count+1
    ptemp[i]=ptot/count
    rhotemp[i]=rhottot/count
    #print i,count,ptemp[i],rhotemp[i]
    acs[i]=cs(vars,ptemp[i],rhotemp[i]);
    
print acs.size
print height.size
plt.plot(height/1.0e6,acs/1000)
#ptemp=numpy.float64(data[1:2048,3])

from math import exp, sqrt, pi





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

#compute lambda0
alam0 = numpy.zeros(asize/4)
alami0 = numpy.zeros(asize/4)
print asize
for i in range(1,(asize/4)-1):
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
    lamdash0[i]=-diff5p(alami0,i,h)
    if lamdash0[i]<-200:
        lamdash0[i]=-200
    count=1
    lamtot=0
    for j in range(i,i+30):
        if j<=(asize/4)-30:
            lamtot=lamtot+lamdash0[j]
            count=count+1
    lamtot=lamtot/count
    lamdash0[i]=lamtot 
    
    
fd = open("cutoff.csv", "w")
for i in range(2,(asize/4)-15):
    h=height[i]-height[i+1]
    #lamdash0[i]=diff3p(alam0,i,h)
    #lamdash0[i]=ptemp[i]/rhotemp[i]
    #lamdash0[i]=diff5p(alami0,i,h)
    #if lamdash0[i]<-200:
    #    lamdash0[i]=-200
    #print h,alam0[i],lamdash0[i]
    #print h,ptemp[i],rhotemp[i],lamdash0[i]
    #atc0[i]=1.0/((vars.fgamma*vars.ggg/(4*pi*acs[i]))*math.sqrt(1+2*lamdash0[i]))
    #atc0[i]=-1.0/((vars.fgamma*vars.ggg/(4*pi*acs[i])))*math.sqrt(1+2*dp_over_rhog(vars, height[i]))
    #atc0[i]=math.sqrt(-1.0/((vars.fgamma*vars.ggg/(2*acs[i]*acs[i])))*(1-vars.fgamma*vars.ggg*dp_over_rhog(vars, height[i])))
    #atc0[i]=math.sqrt(-1.0/((vars.fgamma*vars.ggg/(2*acs[i]*acs[i]))))
    #print i,height[i]/1.0e6,atc0[i]
    #lamdash0[i]=0
    atc0[i]=1.0/((vars.fgamma*vars.ggg/(4*pi*acs[i]))*math.sqrt(1+2*lamdash0[i]))
    #print height[i],atc0[i]
    
    fouts=repr(height[i])+' '+repr(atc0[i])+'\n'
    print(fouts)
    fd.write(fouts)

    

plt.plot(height/1.0e6,atc0)
fd.close()



#plt.plot(height/1.0e6,alam0)
#plt.plot(height/1.0e6,lamdash0)

