EXE=${GRASP}/bin/rtransition_phase
LIBS=-L ${GRASP}/lib/ -ldvd90 -lrang90 -lmcp90 -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/libdvd90 -I ${GRASP}/src/lib/librang90 -I ${GRASP}/src/lib/libmcp90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

OBJS= \
	alclla.o \
	alclla_I.o \
	alcnma.o \
	alcnma_I.o \
	alcnsa.o \
	alcnsa_I.o \
	alcnta.o \
	alcnta_I.o \
	angdata.o \
	angdata_I.o \
	bessj.o \
	bessj_I.o \
	brkt.o \
	brkt_I.o \
	connect.o \
	connect_I.o \
	cpmix.o \
	cpmix_I.o \
	spme.o \
	spme_I.o \
	csfm.o \
	csfm_I.o \
	engout1.o \
	engout1_I.o \
	fname.o \
	fname_I.o \
	lodrwff.o \
	lodrwff_I.o \
	lodrwfi.o \
	lodrwfi_I.o \
	getrmp.o \
	getrmp_I.o \
	getosd.o \
	getosd_I.o \
	iqr.o \
	iqr_I.o \
	isparr.o \
	isparr_I.o \
	itjpor.o \
	itjpor_I.o \
	jcupr.o \
	jcupr_I.o \
	jqsr.o \
	jqsr_I.o \
	ldcsl1.o \
	ldcsl1_I.o \
	ldcsl2.o \
	ldcsl2_I.o \
	ldlbl1.o \
	ldlbl1_I.o \
	ldlbl2.o \
	ldlbl2_I.o \
	lodcslm.o \
	lodcslm_I.o \
	mctin.o \
	mctin_I.o \
	trsort.o \
	trsort_I.o \
	mctout_gg.o \
	mctout_gg_I.o \
	merg12.o \
	merg12_I.o \
	mrgcsl.o \
	mrgcsl_I.o \
	readmix.o \
	readmix_I.o \
	printaLS.o \
	printaLS_I.o \
	printa.o \
	printa_I.o \
	oscl.o \
	oscl_I.o \
	setcslm.o \
	setcslm_I.o \
	strsum.o \
	strsum_I.o \
	testmix.o \
	testmix_I.o \
	bioscl.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
