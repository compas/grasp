.SUFFIXES: .f90 .mod

EXE = rsms
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

APP_OBJ= sms1_C.o teilst_C.o dvpot_C.o\
        sms_I.o      getmixblock_I.o  rintiso_I.o  smsmcp_I.o \
        densmcp_I.o  getsmd_I.o       setdbg_I.o   smsnew_I.o \
        densnew_I.o  polint_I.o       setmcp_I.o   strsum_I.o \
        gco_I.o      rintdens_I.o     setsum_I.o   vinti_I.o \
	engout_I.o wghtd5_I.o \
\
        densmcp.o densnew.o gco.o getsmd.o polint.o rintdens.o \
        rintiso.o setdbg.o setmcp.o setsum.o sms.o sms92.o smsmcp.o \
        smsnew.o strsum.o vinti.o getmixblock.o \
	engout.o wghtd5.o

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
