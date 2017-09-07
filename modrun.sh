#!/bin/bash

# Minimalistic script for run the Unified EMEP model
GRID=EECCA
NLEV=20lev
cd ~/work/EMEP_MSC-W_model.rv4.15.OpenSource/Base_${GRID}_${NLEV}

# Link the input data
inputdir=../
ln -s $inputdir/input/* .        # input files except meteorology
ln -s $inputdir/input/$NLEV/* .  # num-level dependant input
ln -s $inputdir/input/$GRID/* .  # grid dependant input

# Run the model
mpiexec $inputdir/code/Unimod

# Clean the links to the input data
find -maxdepth 1 -type l -delete
