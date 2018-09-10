
.SUFFIXES: .f90 .mod

EXE = rwfnestimate
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ../../lib
MODDIR = ${SRCLIBDIR}/libmod
MODL9290 = ${SRCLIBDIR}/lib9290
MODLRANG90 = ${SRCLIBDIR}/librang90
GRASPLIBS =   -l9290 -lmod  -lrang90

APP_LIBS = -L${GRASPLIB} ${GRASPLIBS} -llapack -lblas

APP_OBJ = \
genrwf_I.o  screenpar_I.o  strsum_I.o  wrtrwf_I.o \
frmhyd_I.o  getinf_I.o  setdbg_I.o  summry_I.o frmrwf_I.o  prtrem_I.o \
setsum_I.o  tail_I.o    frmtfp_I.o  sbstep_I.o  solvh_I.o   tfpot_I.o \
 \
 erwf.o       frmhyd.o     frmrwf.o     frmtfp.o     genrwf.o \
getinfo.o     prtrem.o     sbstep.o     screenpar.o  setdbg.o     setsum.o \
solvh.o      strsum.o     summry.o     tail.o       tfpot.o      wrtrwf.o

$(EXE): ${APP_IOBJ} $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ)  $(APP_LIBS) ${LAPACK_LIBS}

.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I $(MODDIR) -I . -I ${MODL9290} \
	 -I ${MODLRANG90}  -o $@


.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core *.mod

