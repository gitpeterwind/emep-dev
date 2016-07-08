#!/usr/bin/env python
from pylab import *


z=array([ 0.0001, 1.0, 5.0, 10.0 ])
f=linspace(0,10.0,100)  # Fine res
prof1=zeros(len(z))
prof2=zeros(len(f))
z0=0.1

for i in range(0,len(z)):
  prof1[i]= log(z[i]/z0)

for i in range(0,len(f)):
  prof2[i]= log(f[i]/z0)

plot(prof1,z)
plot(prof2,f)
show()
