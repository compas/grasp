EXE=${GRASP}/bin/rangular_mpi
LIBS=-L ${GRASP}/lib/ -lrang90 -lmcp90 -lmpiu90 -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/librang90 -I ${GRASP}/src/lib/libmcp90 -I ${GRASP}/src/lib/mpi90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

OBJS= \
	fndbeg.o \
	fndbeg_I.o \
	getinf.o \
	getinf_I.o \
	outsdampi.o \
	outsdampi_I.o \
	setdbg.o \
	setdbg_I.o \
	setdbgmpi.o \
	setdbgmpi_I.o \
	setmcp.o \
	setmcp_I.o \
	setmcpmpi.o \
	setmcpmpi_I.o \
	setsda.o \
	setsda_I.o \
	setsum.o \
	setsum_I.o \
	settmp.o \
	settmp_I.o \
	sort.o \
	sort_I.o \
	strsum.o \
	strsum_I.o \
	mcpmpi_gg.o \
	mcpmpi_gg_I.o \
	genmcpmpi.o

$(EXE): $(OBJS)
	$(FC_MPI) -o $@ $? $(FC_MPILD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC_MPI) -c $(FC_MPIFLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
