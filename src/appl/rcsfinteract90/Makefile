EXE=${GRASP}/bin/rcsfinteract
LIBS=-L ${GRASP}/lib/ -lrang90 -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/librang90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

OBJS= \
	onescalar1INT.o \
	onescalar1INT_I.o \
	onescalar2INT.o \
	onescalar2INT_I.o \
	onescalarINT.o \
	onescalarINT_I.o \
	el1INT.o \
	el1INT_I.o \
	el2INT.o \
	el2INT_I.o \
	el31INT.o \
	el31INT_I.o \
	el32INT.o \
	el32INT_I.o \
	el33INT.o \
	el33INT_I.o \
	el3INT.o \
	el3INT_I.o \
	el41INT.o \
	el41INT_I.o \
	el4INT.o \
	el4INT_I.o \
	el51INT.o \
	el51INT_I.o \
	el52INT.o \
	el52INT_I.o \
	el53INT.o \
	el53INT_I.o \
	el5INT.o \
	el5INT_I.o \
	Interact_MR.o \
	Interact_MR_I.o \
	Interact_csf.o \
	Interact_csf_I.o \
	getinf.o \
	getinf_I.o \
	lodcsl_CSF.o \
	lodcsl_CSF_I.o \
	lodcsl_MR.o \
	lodcsl_MR_I.o \
	set_CSF_list.o \
	set_CSF_list_I.o \
	set_CSF_number.o \
	set_CSF_number_I.o \
	RCSFinteract.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
