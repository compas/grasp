EXE=${GRASP}/bin/rangular
LIBS=-L ${GRASP}/lib/ -lrang90 -lmcp90 -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/librang90 -I ${GRASP}/src/lib/libmcp90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

OBJS= \
	allocCheck.o \
	allocCheck_I.o \
	fndbeg.o \
	fndbeg_I.o \
	getinf.o \
	getinf_I.o \
	sort.o \
	sort_I.o \
	sortmem.o \
	sortmem_I.o \
	outsda.o \
	outsda_I.o \
	setsda.o \
	setsda_I.o \
	mcp_gg.o \
	mcp_gg_I.o \
	setdbg.o \
	setdbg_I.o \
	setmcp.o \
	setmcp_I.o \
	setmcp2.o \
	setmcp2_I.o \
	setsum.o \
	setsum_I.o \
	settmpGG.o \
	settmpGG_I.o \
	strsum.o \
	strsum_I.o \
	genmcp.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-rm -f $(EXE)
	-rm -f *.o *.mod
