EXE=${GRASP}/bin/rcsfzerofirst
LIBS=-L ${GRASP}/lib/ -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

OBJS= \
	lodcsl_Part.o \
	lodcsl_Part_I.o \
	lodcsl_Zero.o \
	lodcsl_Zero_I.o \
	set_CSF_number.o \
	set_CSF_number_I.o \
	set_CSF_ZFlist.o \
	set_CSF_ZFlist_I.o \
	RCSFzerofirst.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-rm -f $(EXE)
	-rm -f *.o *.mod
