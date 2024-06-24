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
    
nsnap=15

for ((n=0; n <= $nsnap; n++)); 
do

    if [ $n -lt 10 ]; then
	strn="000"$n
    elif [ $n -ge  10 -a $n -lt  100 ]; then
	strn="00"$n
    elif [ $n -ge 100 -a $n -lt 1000 ]; then
	strn="0"$n
    elif [ $n -ge 1000 ]; then
	strn=$n
    else
	echo "too many snapshots."
    fi
    echo $strn

    # Replace the snapshot line in the input file
    sed "s/snapshot=.*/snapshot='$strn'/" input_orig.in > input.in
    if [ ! -d output ]; then
	mkdir output
    fi

    # Run the code with the correct input file for the snapshot
    ./source/flux_feautrier.x
    
    # Rewrite the output snapshot if it already exists
    datadir="output_snapshot"$strn
    if [ -d $datadir ]; then
	rm -rf $datadir
    fi
    mv output $datadir
done



