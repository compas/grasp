.SUFFIXES: .f90 .mod

EXE = rnucleus
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ../../lib
MODDIR = ${SRCLIBDIR}/libmod
MODL9290 = ${SRCLIBDIR}/lib9290
GRASPLIBS =-l9290 -lmod

APP_LIBS = -L${GRASPLIB} ${GRASPLIBS}

APP_OBJ= \
	nucleus_m.o estrms_I.o getcpr_I.o skfun_I.o  \
	estrms.o geniso.o getcpr.o skfun.o

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(iFC_LD) $(APP_OBJ) $(APP_LIBS)

.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I $(MODDIR) -I . -I ${MODL9290} -o $@



clean:
	-rm -f *.o core  *.mod
