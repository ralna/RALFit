#!/bin/bash

if [ $# -eq 0 ]
 then
   file="ARGAUSS"
else
   file=$1    
fi

echo "*************************************"
echo "**                                 **"
echo "**        R A L _ N L L S          **"
echo "**                                 **"
echo "*************************************"

echo " "
echo " ~~~~~~~ $file ~~~~~~~~ "
echo " "
make gen_problem SIFFILE=../cutest/sif/$file
make all
make cutest
./cutest
