#!/usr/bin/env python
import numpy as np
import optparse
import os
from numpy import nan
import matplotlib.mlab  as m
import matplotlib.pylab as plt

# --------- parse inputs -----------------
parser=optparse.OptionParser()
#parser.add_option('-i','--ifile',help="Input file name")

parser.add_option( '-i' ,help="Input file name")
parser.add_option( '-o' ,help="Output file name",default="Screen")
parser.add_option( '-t' ,help="Run on fake test data")

parser.add_option( '--xrange' , help="x-range of data, eg --xrange 50 80", dest='xrange', action='store',nargs=2)
parser.add_option( '--yrange' , help="y-range of data, eg --yrange 50 80", dest='yrange', action='store',nargs=2)

parser.add_option( '--skip' ,help="No. header lines to skip",default=0)
parser.add_option( '-c' ,help="Coord axis")
parser.add_option( '-z', '--z' ,help="Vert axis?",default='z')

opts, args = parser.parse_args()

#opts.t=True
#opts.i="fakedata.txt"
if opts.t:
      print "\n TESTING from script? \n We use fake data.\n"
      ifile=file("fakedata.txt","w")
      ifile.write("Abc z  val\nAAA 1.0  2.3\nBBB 2.0 3.4\nCCC 3.0 5.5\n");
      ifile.close()

if opts.i is None:
      print "\n ERROR No input file specified!!\n"
      parser.print_help()
      exit(-1)
else:
      ifile = opts.i
      print "\n Input file:\n", ifile

#PQR#if opts.xrange is None: opts.xrange = ( -20, 50 )
#PQR#if opts.yrange is None: opts.yrange = (  30, 70 )
#PQR
# --------- Data --------------------------
# recfromtxt: variant of genfromtxt, gets record arrays
# Can access, through e.g. r['z'] or r.z

#Aug2013 r=np.recfromtxt( opts.i ,names=True,dtype=None)
# When names=True, comments are examined anyway for headers.
# From : http://docs.scipy.org/doc/numpy/user/basics.io.genfromtxt.html
# Need to skip lines:

nskip=0
if opts.skip:
  nskip = int(opts.skip)
print "SKIPS ", nskip
#r=np.recfromtxt( opts.i ,names=True,skip_header=nskip,comments="#",dtype=None)
try:
  r=np.recfromtxt( opts.i ,names=True,skip_header=nskip,comments="#",dtype=None)
except:
  print "ERROR reading file. Needs --skip maybe?"
  exit()

print "R=", r
hdrs=r.dtype.names
nrows=len(r)
nhdrs=len(hdrs)

print "recfromtxt: First, last header, size ", r.z, nhdrs, 'x', nrows
print "recfromtxt: hdrs  ", hdrs

#First, we find coord axis

cx = 'z'   # default
if opts.c:
   cx = opts.c   # default

if hdrs.index(cx):
     print "Coord axis found ", cx
     print r[cx]

maxv = -9.99e10   # search for max value
for i in range(0,nhdrs):
   h=hdrs[i]
   #if h == xlab:
   #   xlab=h
   ##elif h == 'x':
   #   xlab='x', x=r[h]
   #elif h == 't':
   #   xlab='t', x=r[h]
   #else:
   if h == cx:
        continue
   if r.dtype[i].kind in ['i','f']:
        print "FLOAT ", h
        if( cx == 'z') :
          plt.plot(r[h], r[cx], label='%s'% hdrs[i] )
        else:
          plt.plot(r[cx], r[h], label='%s'% hdrs[i] )
        if( max( r[h] )  > maxv ) : maxv = max( r[h] )
#   print "for coord axis, we use  ", x

# Dirty plotting trick. The plot looks odd if the max x values is the same as
#  the axis. (Happens for Kz for simple testing). We add a little
v=np.array( plt.axis() )
print "VVVV ", v, "MAX ", maxv
if ( maxv == v[1] ): v[1]=maxv*1.01
plt.axis(v)

if opts.xrange:
   print "xlim", opts.xrange
   x0 = float(opts.xrange[0])
   x1 = float(opts.xrange[1])
   plt.xlim( x0, x1 )

plt.ylabel('z (m)')
plt.legend()
plt.show()


if ifile == "fakedata.txt":
  os.remove(ifile)
