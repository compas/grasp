.SUFFIXES: .f90 .mod

EXE = rangular_mpi
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ../../lib
MODDIR = ${SRCLIBDIR}/libmod
MODL9290 = ${SRCLIBDIR}/lib9290
MODLRANG90 = ${SRCLIBDIR}/librang90
MODLMCP90 = ${SRCLIBDIR}/libmcp90
MODLMPIU90 = ${SRCLIBDIR}/mpi90
GRASPLIBS = -l9290 -lrang90 -l9290 -lmpiu90 -lmcp90 -l9290 -lmod

APP_LIBS = -L ${GRASPLIB} ${GRASPLIBS}


#   Define data types
VASTO = ${MODDIR}/vast_kind_param_M.o

APP_OBJ= \
	fndbeg_I.o mcpmpi_gg_I.o outsdampi_I.o setmcpmpi_I.o setmcp_I.o  settmp_I.o \
	setdbgmpi_I.o  setsda_I.o  sort_I.o \
	setdbg_I.o \
	getinf_I.o setsum_I.o  strsum_I.o \
\
        genmcpmpi.o  mcpmpi_gg.o outsdampi.o setmcpmpi.o \
        fndbeg.o getinf.o setdbgmpi.o setmcp.o setsda.o setsum.o settmp.o \
	setdbg.o \
        sort.o strsum.o

$(EXE): $(APP_OBJ)
	$(FC_MPI) -o $(BINFILE) $(FC_MPILD) $(APP_OBJ) $(APP_LIBS) \
                     ${LAPACK_LIBS}

.f90.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) $< -I $(MODDIR) -I $(MODL9290) -I $(MODLRANG90) -I $(MODLMCP90) -I ${MODDIR} -I $(MODLMPIU90) -o $@


.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core *.mod
