.SUFFIXES: .f90 .mod

EXE = rangular
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ../../lib
MODDIR = ${SRCLIBDIR}/libmod
MODL9290 = ${SRCLIBDIR}/lib9290
MODLRANG90 = ${SRCLIBDIR}/librang90
MODLMCP90 = ${SRCLIBDIR}/libmcp90
GRASPLIBS =  -l9290 -lrang90 -l9290 -lmcp90 -lmod

APP_LIBS = -L ${GRASPLIB} ${GRASPLIBS} 


#   Define data types
VASTO = ${MODDIR}/vast_kind_param_M.o

APP_OBJ= \
	fndbeg_I.o         mcp_gg_I.o      setmcp_I.o  settmpGG_I.o \
	setdbg_I.o	   outsda_I.o   setsda_I.o  sort_I.o \
	getinf_I.o         setmcp2_I.o  setsum_I.o  strsum_I.o \
	allocCheck_I.o sortmem_I.o \
	\
        genmcp.o  mcp_gg.o outsda.o setmcp2.o \
        fndbeg.o getinf.o setdbg.o setmcp.o setsda.o setsum.o settmpGG.o \
        sort.o strsum.o allocCheck.o sortmem.o

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) $(APP_LIBS) $(LAPACK_LIBS)

.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I ${MODDIR} -I $(MODL9290) -I $(MODLRANG90) -I $(MODLMCP90) -I ${MODDIR} -o $@


.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core *.mod

