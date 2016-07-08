#!/bin/bash

# Minimalistic script for run the Unified EMEP model

# Link the input data
inputdir=.
ln -s $inputdir/met/*   .   # Driving meteorology
ln -s $inputdir/input/* .   # Other input files

# Run the model
mpiexec $inputdir/code/Unimod

# Clean the links to the input data
find . -type l -print0 | xargs -0 rm 
