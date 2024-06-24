#!/bin/bash

nsnap=1

if ($nsnap > 10)
   nzeros='000'
if ($nsnap >= 10 and $nsnap < 100 )
   nzeros='00'
if if ($nsnap >= 100 and $nsnap < 1000 )
   nzeros='0'
if ($nsnap >= 1000
   nzeros=''
   
for ((n=0; n <= $nsnap; n++)); 
do    
    sed "s/snapshot=.*/snapshot='$nzeros$n'/" input_orig.in > input.in
    if [ ! -d output ]; then
	mkdir output
    fi
    ./source/flux_feautrier.x
    datadir="output_snapshot"$nzeros$n
    if [ -d $datadir ]; then
	rm -rf $datadir
    fi
    mv output $datadir
done



