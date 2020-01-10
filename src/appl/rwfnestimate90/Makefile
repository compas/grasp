EXE=${GRASP}/bin/rwfnestimate
LIBS=-L ${GRASP}/lib/ -lrang90 -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/librang90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

OBJS= \
	tail.o \
	tail_I.o \
	sbstep.o \
	sbstep_I.o \
	solvh.o \
	solvh_I.o \
	frmhyd.o \
	frmhyd_I.o \
	frmrwf.o \
	frmrwf_I.o \
	frmtfp.o \
	frmtfp_I.o \
	tfpot.o \
	tfpot_I.o \
	prtrem.o \
	prtrem_I.o \
	screenpar.o \
	screenpar_I.o \
	setdbg.o \
	setdbg_I.o \
	setsum.o \
	setsum_I.o \
	strsum.o \
	strsum_I.o \
	summry.o \
	summry_I.o \
	wrtrwf.o \
	wrtrwf_I.o \
	genrwf.o \
	genrwf_I.o \
	getinfo.o \
	getinf_I.o \
	erwf.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-rm -f $(EXE)
	-rm -f *.o *.mod
