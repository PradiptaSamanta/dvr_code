#!/bin/bash

#Set value of maximal order of multipole expansion
lmax=12

#Initialise
make clean all
lprev=0
sed -i 's/.*l_val=.*/  l_val=0/' dvr_diag.f90

#Loop over multipole orders
for ((l=0; l<=lmax; l++)); do
    replacementstring_l="s/l_val=""$lprev""/l_val=""$l""/g"
    perl -pi -e "$replacementstring_l" dvr_diag.f90
    make clean all 
    ./dvr_diag
    wait
    cp transformation_matrix_l""$l"".dat Results/ 
    lprev=$l
done

# Cleanup and reset value of l in code
rm transformation_matrix_l""$lmax"".dat
replacementstring_l="s/para%l=""$lmax""/para%l=0/g"
perl -pi -e "$replacementstring_l" dvr_diag.f90
