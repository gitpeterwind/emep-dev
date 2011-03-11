#!/bin/bash

#Minimum script to run the emep model

#link to the input data
inputdir=.
ln -s $inputdir/input/* .
ln -s $inputdir/met/* .

#define some input data
trendyear=2008 #emission year
runlabel1=Base #short label
runlabel2=Opensource_setup #long label
startdate="2008 01 01" #start date (metdata)
enddate="2008 01 01" #end date (metdata)

#put input data into a temporary file called INPUT.PARA
cat>>    'INPUT.PARA'<<    EOF
$trendyear
$runlabel1
$runlabel2
$startdate
$enddate
EOF

#run the model
mpirun $inputdir/code/Unimod

#clean the links to the input data
ls $inputdir/input|xargs rm
ls $inputdir/met|xargs rm
rm INPUT.PARA
