#!/bin/bash

# Minimalistic script for run the Unified EMEP model

# Link the input data
inputdir=.
ln -s $inputdir/met/*   .   # Driving meteorology
ln -s $inputdir/input/* .   # Other input files

# Define some run parameters
trendyear=2008              # emission year
runlabel1=Base              # short label
runlabel2=Opensource_setup  # long label
startdate="2008 01 01"      # start date (metdata)
  enddate="2008 01 01"      # end date (metdata)

# Put the run parameters in a temporary file
cat > INPUT.PARA << EOF
$trendyear
$runlabel1
$runlabel2
$startdate
$enddate
EOF

# Run the model
mpirun $inputdir/code/Unimod

# Clean the links to the input data and remove INPUT.PARA
find . -type l -print0 | xargs -0 rm 
rm INPUT.PARA
