#!/usr/bin/env python
import sys
import glob

Usage = """

Usage: compinfos.py TestDir

  e.g. compinfos.py ../BM_rv3_7beta19

Intended to be run from one of the BM_ directories that already has an
infov file.  The code will then compare the infosv file in the current
directory with that in the TestDir.
"""

if( len(sys.argv) < 2 ):
  sys.exit("\n\nError!! \n"+Usage)
if( sys.argv[1] == "-h" ):
  sys.exit(Usage)

print "\n\nRunning ", sys.argv[0],  sys.argv[1],"\n"
bfiles=glob.glob('*infov')
tfiles=glob.glob('%s/*infov'% sys.argv[1])

if(len(bfiles)>1):
  sys.exit("ERROR: too many infovs in base dir!")
if(len(tfiles)>1):
  sys.exit("ERROR: too many infovs in test dir!")

bfile=bfiles[0]
tfile=tfiles[0]

def getmaxavg(filename,dmin,davg,dmax):
  try:
    f=file(filename,'r')
  except:
    sys.exit( "ERROR: No input file %s\n"%filename )

  fields=f.readline().split()
  try:
    ivar=fields.index('Name')
  except:
    ivar=fields.index('Parameter')
  imin=fields.index('Minimum')
  iavg=fields.index('Mean')
  imax=fields.index('Maximum')
  for line in f:   
    fields = line.split()        # Temp, we had replaced NO3_f with NO3_F, _c with _C
    key = fields[ivar].upper()   #  and want to still match the values.
    dmin[key] = fields[imin]
    davg[key] = fields[iavg]
    dmax[key] = fields[imax]


(bmin,bavg,bmax) = ({},{},{})
getmaxavg(bfile,bmin,bavg,bmax)

(tmin,tavg,tmax) = ({},{},{})
getmaxavg(tfile,tmin,tavg,tmax)

# We use bmax to obtain the keys, but compare max,min and avg for
# differences
print "-------------------------\n"
Ndiffs=0
Nmissing=0
for key, value in bmax.items():
  try:
    tested = tmax[key] # Just see if we have this in try.

    # if(value != tmax[key] or bavg[key] != tavg[key] or
          # bmin[key] != tmin[key] ):
    if(any((value!=tmax[key],bavg[key]!=tavg[key],bmin[key]!=tmin[key]))):
      print "%-20s Max %s  %s   Avg %s  %s   Min %s  %s"% \
                       (key, 
      bmax[key], tmax[key], 
      bavg[key], tavg[key], 
      bmin[key], tmin[key])
      Ndiffs += 1
  except:
    print "%s doesn't exist in test" % key
    Nmissing += 1

print "\n No. different lines = %d  " % Ndiffs
print "\n No. missing keys    = %d  " % Nmissing
print "-------------------------\n\n"
