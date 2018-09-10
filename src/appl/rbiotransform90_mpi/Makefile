.SUFFIXES: .f90 .mod

EXE = rbiotransform_mpi
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

APP_OBJ=  \
	angdatampi_I.o  getmixmpi_I.o   invmat_I.o   pamtmt_I.o   setcslbmpi_I.o \
	setcslampi_I.o biotr1_I.o   getsmpi_I.o     kapdata_I.o  prsym_I.o    setvec_I.o \
	bndinv_I.o   lodrwffmpi_I.o  qqsortmpi_I.o   tcsl_I.o \
	brkt_I.o     ielsum_I.o   lodrwfimpi_I.o  radfilempi_I.o  tiinigmpi_I.o \
	citragmpi_I.o   ifnmnx_I.o   lulu_I.o     radparmpi_I.o   trpmat_I.o \
	copvec_I.o   inprod_I.o   matml4_I.o   rintff_I.o   ulla_I.o \
	fname_I.o    intrpqf_I.o  mcpoutmpi_gg_I.o   rintii_I.o   vecsum_I.o \
	intrpqi_I.o  orbord_I.o   scalve_I.o   wrtmat_I.o \
        mcpinmpi_I.o ti1tv_I.o lodcslBio_I.o \
	biotrmpi.o fname.o intrpqf.o intrpqi.o kapdata.o lodrwffmpi.o lodrwfimpi.o \
        radparmpi.o rintff.o rintii.o tcsl.o brkt.o getsmpi.o wrtmat.o  \
        biotr1.o radfilempi.o ifnmnx.o ielsum.o copvec.o invmat.o \
        ulla.o trpmat.o matml4.o scalve.o setvec.o pamtmt.o bndinv.o lulu.o \
        inprod.o mcpinmpi.o getmixmpi.o qqsortmpi.o citragmpi.o tiinigmpi.o \
        ti1tv.o vecsum.o prsym.o setcslampi.o setcslbmpi.o angdatampi.o mcpoutmpi_gg.o lodcslBio.o

$(EXE): $(APP_OBJ)
	$(FC_MPI) -o $(BINFILE) $(FC_MPILD) $(APP_OBJ) $(APP_LIBS) \
		$(APP_LIBS)  $(LAPACK_LIBS)

.f90.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) $< -I $(MODDIR) -I . -I ${MODL92} -I ${MODLRANG90} -I $(MODLMCP90)  -I $(MODLMPIU90) -o $@

.f.o:
	$(FC_MPI) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core *.mod


