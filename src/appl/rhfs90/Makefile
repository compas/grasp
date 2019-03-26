.SUFFIXES: .f90 .mod

EXE = rhfs
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ../../lib
MODDIR = ${SRCLIBDIR}/libmod
MODL92 = ${SRCLIBDIR}/lib9290
MODLRANG90 = ${SRCLIBDIR}/librang90
MODLMCP90 = ${SRCLIBDIR}/libmcp90
GRASPLIBS =-l9290 -lmod  -lrang90 -lmcp90

APP_LIBS = -L ${GRASPLIB} ${GRASPLIBS}

APP_OBJ= opt6_C.o\
        getmixblock_I.o  rinthf_I.o  strsum_I.o \
        engouth_I.o  hfsgg_I.o          setdbg_I.o \
        gethfd_I.o   matelt_I.o       setsum_I.o \
 \
        engouth.o gethfd.o hfsgg.o hfs92.o matelt.o  rinthf.o \
        setdbg.o  setsum.o strsum.o getmixblock.o

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) $(APP_LIBS) \
		 $(APP_LIBS) ${LAPACK_LIBS}
.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I $(MODDIR) -I ${MODL92} -I $(MODLRANG90) -I $(MODLMCP90) \
	        -I $(MODDIR) -o $@

.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core  *.mod
