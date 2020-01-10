EXE=${GRASP}/bin/rnucleus
LIBS=-L ${GRASP}/lib/ -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

OBJS= \
	skfun.o \
	skfun_I.o \
	estrms.o \
	estrms_I.o \
	getcpr.o \
	getcpr_I.o \
	geniso.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-rm -f $(EXE)
	-rm -f *.o *.mod
