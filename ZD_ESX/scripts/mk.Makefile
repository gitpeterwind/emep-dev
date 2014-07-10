
# Script to build a makefile for a desired driver (e.g. esx_tester.f90)
# Calls Hugh Pumphrey's fmkmf perl script, using gfortran. Here
# we use the pedantic gfortran compilation. Change as desired, but
# keep the -fdefault-real-8 which provides double precision.

# Test if any argument supplied:
# (-z will test if the expansion of "$1" is a null string or not)
if [ -z "$1" ]; then
  echo "No argument supplied! "
  echo " "
  echo " Usage: mk.Makefile testfile.f90  [MAKEFILE_name]" 
  exit 1
fi

if [ -z "$2" ]; then
  Mfile=Makefile # defaul
  echo "HERE A ",$2
else
  Mfile=$2
  echo "HERE B ",$2
fi

# Now, we need to know where the scripts directory is (home of this)
# and from this, we deduce the directory of ESX_CanChem
#progdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#progdir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
progdir=$( cd "$( dirname "$0" )" && pwd )
topdir=`dirname $progdir`
current="$PWD"

# Use perl or python to solve the tricky problem of getting relative paths
# see stackoverflow.com/questions/2564634/bash-convert-absolute-path-into-relative-path-given-a-current-directory'
# Python comes out cleaner

#reldir=$(perl -MFile::Spec -e 'print File::Spec->abs2rel("'$topdir'","'$current'")')
#reldir2=$(perl -e 'use File::Spec; print File::Spec->abs2rel("'$topdir'","'$current'")')
reldir=$(python -c "import os.path; print os.path.relpath('$topdir','$current')")


echo "SCRIPT DIR ", $progdir 
echo "current DIR ", $current 
echo " -> topdir ", $topdir
echo " -> reldir ", $reldir


# Run perl script, including DO3SE directory and allowing for different fortran suffixes

$progdir/fmkmf  \
     -p .:$reldir/DO3SE/src \
     -f90 "gfortran  -pedantic -Wall -fbounds-check -fdefault-real-8 -finit-real=nan"\
     -tag "(f95|F95|f90|F90)" $1 > $Mfile

echo "(IGNORE warning about iso_fortran_env. This has no impact)"
echo " "
echo "gfortran Makefile created for $1"


