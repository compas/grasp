#!/bin/bash

rnucleus << eof
4.
9.
n
9.
1
1
1
eof

rcsfgenerate << eof
*
0
1s(2,i)2s(2,*)

4s,4p,4d,4f
0,0
2
n
eof

cp rcsf.out rcsf.inp 

rwfnestimate << eof
y
2
*
eof

mv rcsf.inp n4.c
mv rwfn.inp n4.w

rci << eof
y
n4
n
n
n
n
n
n
1
eof

rdensity << eof
y
n4
y
2
eof

cp n4.nw n4NO.w
cp n4.c n4NO.c

rci << eof
y
n4NO
n
n
n
n
n
n
1
eof
echo ' '
echo ' '
echo ' '

echo '------------------------'
echo '   Thomas-Fermi basis'
echo '------------------------'
echo ' '
echo ' '
echo ' '
rmixextract << eof
n4
y
0
n
eof
echo ' '
echo ' '
echo ' '
echo ' '
echo ' '
echo ' '
echo '------------------------'
echo ' Natural orbital basis'
echo '------------------------'

echo ' '
echo ' '
echo ' '
rmixextract << eof
n4NO
y
0
n
eof
