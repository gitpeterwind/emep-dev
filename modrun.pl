#!/bin/bash

#Minimum script to run the emep model

#link the input data
inputdir=./Base
ln -s $inputdir/* .

#define some input data
trendyear=2008 #emission year
runlabel1=rv3_7 #short label
runlabel2=Opensource_setup #long label
startyear=2008 #start year (metdata)
startmonth=1 #start month (metdata)
startday=1 #start day (metdata)
endmonth=1 #end month (metdata)
endday=1 #end day (metdata)

#put input data into a file called INPUT.PARA
cat>>    'INPUT.PARA'<<    EOF
$trendyear
$runlabel1
$runlabel2
$startyear
$startmonth
$startday
$endmonth
$endday
EOF

#run the model
mpirun Unimod

#clean the links to the input data
find . -type l -delete
rm INPUT.PARA
