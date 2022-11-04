# rdensity
Program to compute the density and natural orbitals.


# Installation
Store the rdensity folder in your grasp src/app/ folder, along with all other grasp programs.
Follow the grasp installation procedure.
Look into the Makefile to insert your own set up e.g., 'GRASP' should contain the path to the main grasp directory.
Then type the following commands

>>cd src/app/rdensity
>>make clean
>>make

for compilers options, please look at the Grasp manual.

#Input/output files
Input file:  isodata, name.c, name.w, name.(c)m
Output files: name.nw, name.(c)d

the file name.nw contains the natural orbitals and the file name.(c)d contains the radial density.
