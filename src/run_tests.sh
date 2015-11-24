#!/bin/bash

#if [ $# -eq 0 ]
# then
#   FNAME="results.out"
#else
#   FNAME=$1    
#fi

echo "*************************************"
echo "**                                 **"
echo "**        R A L _ N L L S          **"
echo "**                                 **"
echo "*************************************"

> results.out
> results.out_data
cat sif_names.txt | while read file; do
    echo " "
    echo " ~~~~~~~ $file ~~~~~~~~ "
    echo " "
    make gen_problem SIFFILE=../cutest/sif/$file
    make all
    make cutest
    ./cutest
done

if [ $# -eq 1 ]
 then
    mv results.out $1
    mv results.out_data $1_data
fi
