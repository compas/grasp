.SUFFIXES: .f90 .mod

EXE = rtransition_mpi
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ../../lib
MODDIR = ${SRCLIBDIR}/libmod
MODL92 = ${SRCLIBDIR}/lib9290
MODLRANG90 = ${SRCLIBDIR}/librang90
MODLMCP90 = ${SRCLIBDIR}/libmcp90
MODLMPIU90 = ${SRCLIBDIR}/mpi90

GRASPLIBS =-l9290 -lmod -lrang90 -lmcp90 -ldvd90 -lmpiu90

APP_LIBS = -L${GRASPLIB} ${GRASPLIBS}

APP_OBJ= \
	fname_I.o ldcsl1mpi_I.o ldcsl2mpi_I.o lodcslm_I.o  \
	merg12mpi_I.o mrgcslmpi_I.o ldlbl1_I.o ldlbl2_I.o printaLS_I.o \
	readmixmpi_I.o iqr_I.o isparr_I.o itjpor_I.o jcupr_I.o jqsr_I.o lodrwfimpi_I.o \
	lodrwffmpi_I.o getosdmpi_I.o brkt_I.o getrmpmpi_I.o strsum_I.o osclmpi_I.o \
	connect_I.o cpmix_I.o mctoutmpi_gg_I.o\
	alcnsa_I.o alcnta_I.o  mctinmpi_I.o csfm_I.o printa_I.o spme_I.o trsortmpi_I.o \
	bessj_I.o alclla_I.o alcnma_I.o setcslm_I.o angdatampi_I.o engout1_I.o testmix_I.o \
	\
	biosclmpi.o fname.o ldcsl1mpi.o ldcsl2mpi.o lodcslm.o  \
	merg12mpi.o mrgcslmpi.o ldlbl1.o ldlbl2.o printaLS.o \
	readmixmpi.o iqr.o isparr.o itjpor.o jcupr.o jqsr.o lodrwfimpi.o \
	lodrwffmpi.o getosdmpi.o brkt.o getrmpmpi.o strsum.o osclmpi.o \
	connect.o cpmix.o mctoutmpi_gg.o\
	alcnsa.o alcnta.o  mctinmpi.o csfm.o printa.o spme.o trsortmpi.o \
	bessj.o alclla.o alcnma.o setcslm.o angdatampi.o engout1.o testmix.o

$(EXE): $(APP_OBJ)
	$(FC_MPI) -o $(BINFILE) $(FC_MPILD) $(APP_OBJ) $(APP_LIBS) \
		$(APP_LIBS) ${LAPACK_LIBS} 

.f90.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) $< -I $(MODDIR) -I . -I ${MODL92} -I ${MODLRANG90} -I $(MODLMCP90) -I $(MODLMPIU90) -o $@

.f.o:
	$(FC_MPI) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core *.mod

