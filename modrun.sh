#!/bin/bash
# Basic script to run the emep model

# Link the input data into the run directory
inputdir=./input_data
ln -s $inputdir/* .

# Define the run parameters
trendyear=2008              # emission year
runlabel1=rv3_7             # short label
runlabel2=Opensource_setup  # long label
startdate="2008 01 01"      # start date (metdata)
  enddate="2008 12 31"      # end date (metdata)

# Put run parameters into a temporary file
echo -e "$trendyear\n$runlabel1\n$runlabel2\n$startdate\n$enddate" > INPUT.PARA

# Run the model
mpirun Unimod

# Clean the links to the input data
for L in *; do
  test -L $L && rm -f $L
done
rm INPUT.PARA
