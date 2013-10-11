#!/usr/bin/env python
from pylab import *
import optparse

# --------- parse arguments ------------------
parser=optparse.OptionParser()
parser.add_option('-i',help="input file name")
parser.add_option('-c',help="Component name, e.g. O3")
parser.add_option('-o',help="Output file name (e.g. plotO3.png, plotOH.eps), otherwise print to screen")
parser.add_option('-u',help="Undefined, eg -u 999 or -u 1.0e90  - only plot values below this")
parser.add_option('--tfmt',help="Format string for time labels, e.g. t=%8.3f", dest='tfmt',default='%8.3f')
parser.add_option('--xlow',help="Min x value")
parser.add_option('--zmax',help="Max z value")
opts, args = parser.parse_args()

if opts.c is None:
   print "Need to specify component or header"
   exit()
comp=opts.c

if opts.i:
  v=loadtxt( opts.i )
else:
  comp="OH"
  v=loadtxt( "Zarray_%s.txt" % comp )


# Assumes data has z values in first column and time in first row.
z=v[1:,0]
times=v[0,1:]

Undef= 1.5*v[1:,1:].max() # allows all values
if opts.u: Undef= float( opts.u )
print "plotUNDEF ", Undef
#Test cycling
colors=('k','r','m','c','b','g','r','#aaaaaa')
linestyles=('-','--','-.',':')

figure()
n = 0
for n in range( n, len(times)):

 #c=v[1+n,1:]    # concentrations in z, 1+ needed skips axis
 c=v[1:,1+n]    # concentrations in z, 1+ needed skips axis
 f=c<Undef  # filter out undefined
 #print "plotCCC n ", n, c

 #print "CMINMAX ", n, min(c), max(c)
 #print "C", c
 #print "U<",Undef,":", c[f]

 #plot(c[f],z[f], label="t=%8.3f" %  times[n] )
 nn=min(n,3) # Keeps color+linestyle index in range
 plot(c[f],z[f], color=colors[nn],ls=linestyles[nn], lw=2,label= opts.tfmt %  times[n] )
 title(comp)
 #axis([0,5,0,10.0])
 #print times[n]
 #print c

[xmin, xmax, ymin, ymax] = axis()
#print "plotAXIS ", axis()
#print "plotDATA ", v[1:,1]
if( opts.zmax ): ymax = float( opts.zmax )
axis([xmin, xmax, ymin, ymax])
# Adjust limits
if opts.xlow: xmin= float( opts.xlow )
xmax = 1.05*xmax 
axis([xmin, xmax, ymin, ymax])
xlabel("Value")
ylabel("Height (m)" )

legend(loc=2)
#plot(t20,z)


if opts.o:
  savefig("%s" % opts.o)
else:
  show()

#plot(x)
