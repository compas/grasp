SCRIPTS="lscomp.pl rsave"
PROGRAMS="
rasfsplit
rcsfblock
rcsfmr
rcsfsplit
rhfs_lsj
rlevelseV
rlevels
rmixaccumulate
rmixextract
rseqenergy
rseqhfs
rseqtrans
rtabhfs
rtablevels
rtabtrans1
rtabtrans2
rtabtransE1
rwfnmchfmcdf
rwfnplot
rwfnrelabel
rwfnrotate
wfnplot
"
# rcsfratip was not being compiled in the original Makefile for some reason.

# Generate Makefile
BINARIES_FLAT="$(for p in $SCRIPTS $PROGRAMS; do echo -n " \${GRASP}/bin/$p"; done)"
cat <<-EOF | sed 's/    /\t/' > Makefile
	LIBS=-L ${GRASP}/lib/ -l9290 -lmod
	FC_MODULES=-I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

	all: ${BINARIES_FLAT}

EOF
for script in ${SCRIPTS}; do
	cat <<-EOF | sed 's/    /\t/' >> Makefile
		\${GRASP}/bin/${script}: ${script}
		    cp \$^ \$@
		    chmod u+x \$@

	EOF
done
for program in ${PROGRAMS}; do
	cat <<-EOF | sed 's/    /\t/' >> Makefile
		\${GRASP}/bin/${program}: ${program}.o
		    \$(FC) -o \$@ \$? \$(FC_LD) \$(LIBS) \$(LAPACK_LIBS)

	EOF
done
cat <<-EOF | sed 's/    /\t/' >> Makefile
	%.o: %.f90
	    \$(FC) -c \$(FC_FLAGS) \$(FC_MODULES) -o \$@ \$<

	clean:
	    -rm -f ${BINARIES_FLAT}
	    -rm -f *.o *.mod
EOF

# Generate CMakeLists
echo -n > CMakeLists.txt
for script in ${SCRIPTS}; do
	cat <<-EOF | sed 's/    /\t/' >> CMakeLists.txt
		install(PROGRAMS ${script} DESTINATION bin/)
	EOF
done
echo >> CMakeLists.txt
for program in ${PROGRAMS}; do
	cat <<-EOF | sed 's/    /\t/' >> CMakeLists.txt
	add_executable(${program} ${program}.f90)
	target_link_libraries_Fortran(${program} PRIVATE mod 9290)
	install(TARGETS ${program} DESTINATION bin/)

	EOF
done
