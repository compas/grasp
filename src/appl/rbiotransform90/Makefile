EXE=${GRASP}/bin/rbiotransform
LIBS=-L ${GRASP}/lib/ -lrang90 -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/librang90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

OBJS= \
	angdata.o \
	angdata_I.o \
	bndinv.o \
	bndinv_I.o \
	brkt.o \
	brkt_I.o \
	wrtmat.o \
	wrtmat_I.o \
	scalve.o \
	scalve_I.o \
	copvec.o \
	copvec_I.o \
	fname.o \
	fname_I.o \
	getmix.o \
	getmix_I.o \
	gets.o \
	gets_I.o \
	ielsum.o \
	ielsum_I.o \
	ifnmnx.o \
	ifnmnx_I.o \
	rintff.o \
	rintff_I.o \
	rintii.o \
	rintii_I.o \
	inprod.o \
	inprod_I.o \
	intrpqf.o \
	intrpqf_I.o \
	intrpqi.o \
	intrpqi_I.o \
	invmat.o \
	invmat_I.o \
	kapdata.o \
	kapdata_I.o \
	lodcslBio.o \
	lodcslBio_I.o \
	lodrwff.o \
	lodrwff_I.o \
	lodrwfi.o \
	lodrwfi_I.o \
	prsym.o \
	prsym_I.o \
	lulu.o \
	lulu_I.o \
	setvec.o \
	setvec_I.o \
	matml4.o \
	matml4_I.o \
	qqsort.o \
	qqsort_I.o \
	mcpin.o \
	mcpin_I.o \
	mcpout_gg.o \
	mcpout_gg_I.o \
	pamtmt.o \
	pamtmt_I.o \
	radfile.o \
	radfile_I.o \
	radpar.o \
	radpar_I.o \
	setcslb.o \
	setcslb_I.o \
	tcsl.o \
	tcsl_I.o \
	ti1tv.o \
	ti1tv_I.o \
	tiinig.o \
	tiinig_I.o \
	trpmat.o \
	trpmat_I.o \
	ulla.o \
	ulla_I.o \
	vecsum.o \
	vecsum_I.o \
	genmcp.o \
	genmcp_I.o \
	citrag.o \
	citrag_I.o \
	biotr1.o \
	biotr1_I.o \
	biotr.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
