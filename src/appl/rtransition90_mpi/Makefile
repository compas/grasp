EXE=${GRASP}/bin/rtransition_mpi
LIBS=-L ${GRASP}/lib/ -lmcp90 -ldvd90 -lrang90 -lmpiu90 -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/libmcp90 -I ${GRASP}/src/lib/libdvd90 -I ${GRASP}/src/lib/librang90 -I ${GRASP}/src/lib/mpi90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

OBJS= \
	alclla.o \
	alclla_I.o \
	alcnma.o \
	alcnma_I.o \
	alcnsa.o \
	alcnsa_I.o \
	alcnta.o \
	alcnta_I.o \
	angdatampi.o \
	angdatampi_I.o \
	bessj.o \
	bessj_I.o \
	brkt.o \
	brkt_I.o \
	connect.o \
	connect_I.o \
	cpmix.o \
	cpmix_I.o \
	spme.o \
	spme_I.o \
	csfm.o \
	csfm_I.o \
	engout1.o \
	engout1_I.o \
	fname.o \
	fname_I.o \
	getrmpmpi.o \
	getrmpmpi_I.o \
	lodrwffmpi.o \
	lodrwffmpi_I.o \
	lodrwfimpi.o \
	lodrwfimpi_I.o \
	getosdmpi.o \
	getosdmpi_I.o \
	iqr.o \
	iqr_I.o \
	isparr.o \
	isparr_I.o \
	itjpor.o \
	itjpor_I.o \
	jcupr.o \
	jcupr_I.o \
	jqsr.o \
	jqsr_I.o \
	ldcsl1mpi.o \
	ldcsl1mpi_I.o \
	ldcsl2mpi.o \
	ldcsl2mpi_I.o \
	ldlbl1.o \
	ldlbl1_I.o \
	ldlbl2.o \
	ldlbl2_I.o \
	lodcslm.o \
	lodcslm_I.o \
	trsortmpi.o \
	trsortmpi_I.o \
	mctinmpi.o \
	mctinmpi_I.o \
	mctoutmpi_gg.o \
	mctoutmpi_gg_I.o \
	merg12mpi.o \
	merg12mpi_I.o \
	mrgcslmpi.o \
	mrgcslmpi_I.o \
	readmixmpi.o \
	readmixmpi_I.o \
	printaLS.o \
	printaLS_I.o \
	printa.o \
	printa_I.o \
	osclmpi.o \
	osclmpi_I.o \
	setcslm.o \
	setcslm_I.o \
	strsum.o \
	strsum_I.o \
	testmix.o \
	testmix_I.o \
	biosclmpi.o

$(EXE): $(OBJS)
	$(FC_MPI) -o $@ $? $(FC_MPILD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC_MPI) -c $(FC_MPIFLAGS) $(FC_MODULES) -o $@ $<

clean:
	-rm -f $(EXE)
	-rm -f *.o *.mod
