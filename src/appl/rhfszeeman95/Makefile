EXE=${GRASP}/bin/rhfszeeman95
LIBS=-L ${GRASP}/lib/ -l9290 -lmod -lrang90 -lmcp90 -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod -I ${GRASP}/src/lib/librang90 -I ${GRASP}/src/lib/libmcp90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

OBJS= \
	opt6_C.o \
	engouth_I.o \
	gethfd_I.o \
	hfszeeman_I.o \
	matelt_I.o \
	rinthf_I.o \
	setdbg_I.o \
	setsum_I.o \
	strsum_I.o \
	getmixblock_I.o \
	engouth.o \
	gethfd.o \
	hfszeeman.o \
	hfszeeman06.o \
	matelt.o \
	rinthf.o \
	setdbg.o \
	setsum.o \
	strsum.o \
	getmixblock.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
