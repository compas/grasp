EXE=${GRASP}/bin/rmcdhf_mem_mpi
LIBS=-L ${GRASP}/lib/ -ldvd90 -lmpiu90 -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/libdvd90 -I ${GRASP}/src/lib/mpi90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

OBJS= \
	rmcdhf_mem_C.o \
	read02_mem.o \
	read02_mem_I.o \
	read03_mem.o \
	read03_mem_I.o \
	read1_mem.o \
	read1_mem_I.o \
	read2_mem.o \
	read2_mem_I.o \
	read3_mem.o \
	read3_mem_I.o \
	systemmem.o \
	systemmem_I.o \
	systemfreemem.o \
	systemfreemem_I.o \
	setmcp_mem.o \
	setmcp_mem_I.o \
	dsubrs.o \
	dsubrs_I.o \
	fco.o \
	fco_I.o \
	gco.o \
	gco_I.o \
	setcof.o \
	setcof_I.o \
	xpot.o \
	xpot_I.o \
	ypot.o \
	ypot_I.o \
	lagcon.o \
	lagcon_I.o \
	dacon.o \
	dacon_I.o \
	cofpotmpi.o \
	cofpotmpi_I.o \
	consis.o \
	consis_I.o \
	csfwgt.o \
	csfwgt_I.o \
	dampck.o \
	dampck_I.o \
	dampor.o \
	dampor_I.o \
	defcor.o \
	defcor_I.o \
	eigen.o \
	eigen_I.o \
	engoutgg.o \
	engoutgg_I.o \
	endsum.o \
	endsum_I.o \
	estim.o \
	estim_I.o \
	getaldwt.o \
	getaldwt_I.o \
	getaldmpi.o \
	getaldmpi_I.o \
	prtrsl.o \
	prtrsl_I.o \
	getoldwt.o \
	getoldwt_I.o \
	getoldmpi.o \
	getoldmpi_I.o \
	getscdmpi.o \
	getscdmpi_I.o \
	hmoutmpi.o \
	hmoutmpi_I.o \
	setxuv.o \
	setxuv_I.o \
	setxv.o \
	setxv_I.o \
	setxz.o \
	setxz_I.o \
	out.o \
	out_I.o \
	in.o \
	in_I.o \
	prwf.o \
	prwf_I.o \
	outbnd.o \
	outbnd_I.o \
	newe.o \
	newe_I.o \
	solve.o \
	solve_I.o \
	setlagmpi.o \
	setlagmpi_I.o \
	orthor.o \
	orthor_I.o \
	orthy.o \
	orthy_I.o \
	improvmpi.o \
	improvmpi_I.o \
	ispar.o \
	ispar_I.o \
	itjpo.o \
	itjpo_I.o \
	lodcslmpiGG.o \
	lodcslmpiGG_I.o \
	lodcsh2GG.o \
	lodcsh2GG_I.o \
	maxarr.o \
	maxarr_I.o \
	newcompi.o \
	newcompi_I.o \
	orbout.o \
	orbout_I.o \
	setcslmpi.o \
	setcslmpi_I.o \
	setdbg.o \
	setdbg_I.o \
	setdbgmpi.o \
	setdbgmpi_I.o \
	setham.o \
	setham_I.o \
	setmcp.o \
	setmcp_I.o \
	setmix.o \
	setmix_I.o \
	setsum.o \
	setsum_I.o \
	strsum.o \
	strsum_I.o \
	maneigmpi.o \
	maneigmpi_I.o \
	matrixmpi.o \
	matrixmpi_I.o \
	scfmpi.o \
	scfmpi_I.o \
	rscfmpivu.o

$(EXE): $(OBJS)
	$(FC_MPI) -o $@ $? $(FC_MPILD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC_MPI) -c $(FC_MPIFLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
