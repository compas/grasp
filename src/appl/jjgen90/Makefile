EXE=${GRASP}/bin/jjgen
LIBS=-L ${GRASP}/lib/ -lmod -l9290 -lrang90
FC_MODULES= -I ${GRASP}/src/lib/libmod -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/librang90

OBJS= \
	lockad.o \
	lockad_I.o \
	reada.o \
	reada_I.o \
	lasa1.o \
	lasa1_I.o \
	lasa2.o \
	lasa2_I.o \
	adder.o \
	adder_I.o \
	slug.o \
	slug_I.o \
	kopp1.o \
	kopp1_I.o \
	kopp2.o \
	kopp2_I.o \
	genb.o \
	gen_I.o \
	sluggo.o \
	sluggo_I.o \
	test.o \
	test_I.o \
	mergeb.o \
	mergeb_I.o \
	blanda.o \
	blanda_I.o \
	blandb.o \
	blandb_I.o \
	blandc.o \
	blandc_I.o \
	copy7t9.o \
	copy7t9_I.o \
	fivefirst.o \
	fivefirst_I.o \
	fivelines.o \
	fivelines_I.o \
	lika.o \
	lika_I.o \
	matain.o \
	matain_I.o \
	matbin.o \
	matbin_I.o \
	matcin.o \
	matcin_I.o \
	merge.o \
	merge_I.o \
	open79.o \
	open79_I.o \
	reffa.o \
	reffa_I.o \
	jjgen15b.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-rm -f $(EXE)
	-rm -f *.o *.mod
