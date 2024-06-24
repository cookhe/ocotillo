#!/bin/bash


# Check if the executable exists
if [ ! -f ./source/flux_feautrier.x ]; then 
    echo "You need to compile and run the code"
    echo "before running this script."
fi

# Save the original input so it is not rewritten
if [ ! -f input_orig.in ]; then
    cp input.in input_orig.in
fi
    
nsnap=1

for ((n=0; n <= $nsnap; n++)); 
do

    if [ $n > 10 ]; then
	nzeros='000'
    elif [ $n >= 10 && $n < 100 ]; then
	nzeros='00'
    elif [ $n >= 100 && $n < 1000 ]; then
	nzeros='0'
    elif [ $n >= 1000 ]; then    
	nzeros=''
    else
	echo "too many snapshots."
    fi

    # Replace the snapshot line in the input file
    sed "s/snapshot=.*/snapshot='$nzeros$n'/" input_orig.in > input.in
    if [ ! -d output ]; then
	mkdir output
    fi

    ./source/flux_feautrier.x
    
# Rewrite the output snapshot if it already exists
    datadir="output_snapshot"$nzeros$n
    if [ -d $datadir ]; then
	rm -rf $datadir
    fi
    mv output $datadir
done



