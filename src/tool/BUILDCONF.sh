SCRIPTS="
rsave
lscomp.pl
rwfnpyplot
"

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
rwfntotxt
fical
"
# rcsfratip was not being compiled in the original ${MAKEFILE} for some reason.

# Generate ${MAKEFILE}
BINARIES_FLAT="$(for p in $SCRIPTS $PROGRAMS; do echo -n " \${GRASP}/bin/$p"; done)"
cat <<-EOF | sed 's/    /\t/' > ${MAKEFILE}
	LIBS=-L \${GRASP}/lib/ -l9290 -lmod
	FC_MODULES=-I \${GRASP}/src/lib/lib9290 -I \${GRASP}/src/lib/libmod

	all: ${BINARIES_FLAT}

EOF
for script in ${SCRIPTS}; do
	cat <<-EOF | sed 's/    /\t/' >> ${MAKEFILE}
		\${GRASP}/bin/${script}: ${script}
		    cp \$^ \$@
		    chmod u+x \$@

	EOF
done
for program in ${PROGRAMS}; do
	cat <<-EOF | sed 's/    /\t/' >> ${MAKEFILE}
		\${GRASP}/bin/${program}: ${program}.o
		    \$(FC) -o \$@ \$? \$(FC_LD) \$(LIBS) \$(LAPACK_LIBS)

	EOF
done
cat <<-EOF | sed 's/    /\t/' >> ${MAKEFILE}
	%.o: %.f90
	    \$(FC) -c \$(FC_FLAGS) \$(FC_MODULES) -o \$@ \$<

	clean:
	    -rm -f ${BINARIES_FLAT}
	    -rm -f *.o *.mod
EOF

# Generate CMakeLists
echo -n > ${CMAKELISTSTXT}
for script in ${SCRIPTS}; do
	cat <<-EOF | sed 's/    /\t/' >> ${CMAKELISTSTXT}
		install(PROGRAMS ${script} DESTINATION bin/)
	EOF
done
echo >> ${CMAKELISTSTXT}
for program in ${PROGRAMS}; do
	cat <<-EOF | sed 's/    /\t/' >> ${CMAKELISTSTXT}
	add_executable(${program} ${program}.f90)
	target_link_libraries_Fortran(${program} PRIVATE mod 9290)
	install(TARGETS ${program} DESTINATION bin/)

	EOF
done
