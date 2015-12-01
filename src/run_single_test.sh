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
cd ../cutest/sif/
runcutest -p ral_nlls --decode $file
echo '$file'_progress.out
mv 'RAL_NLLS_iter.sum' '../../src/data/'$file'_progress.out'
cd ../../src/
