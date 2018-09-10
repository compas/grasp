.SUFFIXES: .f90 .mod

EXE = rcsfinteract
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ../../lib
MODDIR = ${SRCLIBDIR}/libmod
MODL9290 = ${SRCLIBDIR}/lib9290
MODLRANG90 = ${SRCLIBDIR}/librang90
GRASPLIBS =-l9290 -lmod -lrang90

APP_LIBS = -L ${GRASPLIB} ${GRASPLIBS}  

#   Define data types
VASTO = ${MODDIR}/vast_kind_param_M.o

APP_OBJ= \
	getinf_I.o Interact_MR_I.o \
	lodcsl_MR_I.o lodcsl_CSF_I.o \
	set_CSF_list_I.o set_CSF_number_I.o \
	onescalarINT_I.o onescalar1INT_I.o onescalar2INT_I.o \
	el1INT_I.o el2INT_I.o el3INT_I.o el31INT_I.o el32INT_I.o \
	el33INT_I.o el4INT_I.o el41INT_I.o \
	el5INT_I.o el51INT_I.o el52INT_I.o el53INT_I.o Interact_csf_I.o \
\
        RCSFinteract.o getinf.o Interact_MR.o \
        lodcsl_MR.o lodcsl_CSF.o \
	set_CSF_list.o set_CSF_number.o \
	onescalarINT.o onescalar1INT.o onescalar2INT.o \
	el1INT.o el2INT.o el3INT.o el31INT.o el32INT.o el33INT.o \
	el4INT.o el41INT.o \
	el5INT.o el51INT.o el52INT.o el53INT.o Interact_csf.o 

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(APP_OBJ) $(APP_LIBS) \
		$(APP_LIBS) 



.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I $(MODDIR) -I ${MODL9290} -I $(MODLRANG90) -I $(MODDIR) -o $@

.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core  *.mod
