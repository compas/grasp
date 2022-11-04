#!/bin/bash

rnucleus << eof
118.
294.
n
294.
1.
1.
1.
eof

rcsfgenerate << eof
 * ! Orbital order
           6  ! Selected core
5f(14,c)6d(10,c)7s(2,i)7p(6,i)
*
7s,7p,6d,5f
           0           0  ! Lower and higher 2*J
           0  ! Number of excitations
n
eof

mv rcsf.out rcsf.inp

rangular << eof
y
eof

rwfnestimate << eof
y
2
*
eof

rmcdhf << eof
y
1
*
*
100
eof

rsave DHF_Og

rci << eof
y   
DHF_Og
y  
y   
1.E-6
y 
n 
n 
y 
5 
1
eof

rdensity << eof
y
DHF_Og
y
2
eof

rm mcp*
rm rci.res
rm rcsf*
rm clist.new
