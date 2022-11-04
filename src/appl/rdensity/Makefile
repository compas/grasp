EXE=${GRASP}/bin/rdensity
LIBS=-L ${GRASP}/lib/ -l9290 -lmod -lrang90 -lmcp90 -lmod -l9290
FC_MODULES= -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod -I ${GRASP}/src/lib/librang90 -I ${GRASP}/src/lib/libmcp90 -I ${GRASP}/src/lib/libmod -I ${GRASP}/src/lib/lib9290

OBJS= \
	rdensity_C.o \
	teilst_C.o \
	cxk_I.o \
	getmixblock_I.o \
	getsmd_I.o \
	polint_I.o \
	rdensity_cal_I.o \
	rintdens_I.o \
	rintdensvec_I.o \
	natorbnew_I.o \
	setdbg_I.o \
	setdens_I.o \
	rdensity.o \
	cxk.o \
	getmixblock.o \
	getsmd.o \
	polint.o \
	rdensity_cal.o \
	rintdens.o \
	rintdensvec.o \
	natorbnew.o \
	setdbg.o \
	setdens.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
