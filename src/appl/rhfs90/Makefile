EXE=${GRASP}/bin/rhfs
LIBS=-L ${GRASP}/lib/ -lrang90 -lmcp90 -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/librang90 -I ${GRASP}/src/lib/libmcp90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

OBJS= \
	engouth.o \
	engouth_I.o \
	gethfd.o \
	gethfd_I.o \
	getmixblock.o \
	getmixblock_I.o \
	opt6_C.o \
	matelt.o \
	matelt_I.o \
	rinthf.o \
	rinthf_I.o \
	setdbg.o \
	setdbg_I.o \
	setsum.o \
	setsum_I.o \
	strsum.o \
	strsum_I.o \
	hfsgg.o \
	hfsgg_I.o \
	hfs92.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
