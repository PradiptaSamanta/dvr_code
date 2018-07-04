#!/bin/bash

#Set value of maximal order of multipole expansion
lmax=10

#Initialise
make clean all
lprev=0

#Loop over multipole orders
for ((l=0; l<=lmax; l++)); do
    replacementstring_l="s/l_val=""$lprev""/l_val=""$l""/g"
    perl -pi -e "$replacementstring_l" radial_elements.f90
    make clean all 
    ./rad_elements
    wait
    cp rad_elements_l""$l"".dat Results/ 
    lprev=$l
done

# Cleanup and reset value of l in code
rm rad_elements_l""$lmax"".dat
replacementstring_l="s/para%l=""$lmax""/para%l=0/g"
perl -pi -e "$replacementstring_l" radial_elements.f90
