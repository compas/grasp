EXE=${GRASP}/bin/rbiotransform_mpi
LIBS=-L ${GRASP}/lib/ -lrang90 -lmpiu90 -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/librang90 -I ${GRASP}/src/lib/mpi90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

OBJS= \
	angdatampi.o \
	angdatampi_I.o \
	bndinv.o \
	bndinv_I.o \
	brkt.o \
	brkt_I.o \
	vecsum.o \
	vecsum_I.o \
	wrtmat.o \
	wrtmat_I.o \
	copvec.o \
	copvec_I.o \
	fname.o \
	fname_I.o \
	getmixmpi.o \
	getmixmpi_I.o \
	getsmpi.o \
	getsmpi_I.o \
	ielsum.o \
	ielsum_I.o \
	ifnmnx.o \
	ifnmnx_I.o \
	inprod.o \
	inprod_I.o \
	invmat.o \
	invmat_I.o \
	rintff.o \
	rintff_I.o \
	rintii.o \
	rintii_I.o \
	intrpqf.o \
	intrpqf_I.o \
	intrpqi.o \
	intrpqi_I.o \
	kapdata.o \
	kapdata_I.o \
	lodcslBio.o \
	lodcslBio_I.o \
	lodrwffmpi.o \
	lodrwffmpi_I.o \
	lodrwfimpi.o \
	lodrwfimpi_I.o \
	prsym.o \
	prsym_I.o \
	lulu.o \
	lulu_I.o \
	setvec.o \
	setvec_I.o \
	matml4.o \
	matml4_I.o \
	qqsortmpi.o \
	qqsortmpi_I.o \
	mcpinmpi.o \
	mcpinmpi_I.o \
	mcpoutmpi_gg.o \
	mcpoutmpi_gg_I.o \
	pamtmt.o \
	pamtmt_I.o \
	radfilempi.o \
	radfilempi_I.o \
	radparmpi.o \
	radparmpi_I.o \
	scalve.o \
	scalve_I.o \
	setcslampi.o \
	setcslampi_I.o \
	setcslbmpi.o \
	setcslbmpi_I.o \
	tcsl.o \
	tcsl_I.o \
	ti1tv.o \
	ti1tv_I.o \
	tiinigmpi.o \
	tiinigmpi_I.o \
	trpmat.o \
	trpmat_I.o \
	ulla.o \
	ulla_I.o \
	citragmpi.o \
	citragmpi_I.o \
	biotr1.o \
	biotr1_I.o \
	biotrmpi.o

$(EXE): $(OBJS)
	$(FC_MPI) -o $@ $? $(FC_MPILD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC_MPI) -c $(FC_MPIFLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
