.SUFFIXES: .f90 .mod

EXE = rcsfzerofirst
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ../../lib
MODDIR = ${SRCLIBDIR}/libmod
MODL9290 = ${SRCLIBDIR}/lib9290
GRASPLIBS =-l9290 -lmod

APP_LIBS = -L${GRASPLIB} ${GRASPLIBS}

#   Define data types
#VASTO = ${MODDIR}/vast_kind_param_M.o

APP_OBJ= \
	lodcsl_Zero_I.o lodcsl_Part_I.o \
	set_CSF_ZFlist_I.o set_CSF_number_I.o \
\
        RCSFzerofirst.o \
        lodcsl_Zero.o lodcsl_Part.o \
	set_CSF_ZFlist.o set_CSF_number.o

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD)  $(APP_OBJ) $(APP_LIBS) \
 		$(APP_LIBS)

.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I $(MODDIR) -I ${MODL9290} -I $(MODDIR) -o $@

.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core  *.mod
