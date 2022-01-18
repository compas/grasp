EXE=${GRASP}/bin/rmcdhf_mem_mpi
LIBS=-L ${GRASP}/lib/ -ldvd90 -lmpiu90 -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/libdvd90 -I ${GRASP}/src/lib/mpi90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

OBJS = rmcdhf_mem_C.o \
     consis_I.o csfwgt_I.o dacon_I.o dampck_I.o dampor_I.o defcor_I.o \
     dsubrs_I.o eigen_I.o endsum_I.o engoutgg_I.o estim_I.o fco_I.o gco_I.o getaldwt_I.o \
     getoldwt_I.o in_I.o lagcon_I.o maxarr_I.o \
     newe_I.o orbout_I.o \
     orthy_I.o \
     out_I.o outbnd_I.o prtrsl_I.o prwf_I.o setcof_I.o \
     setham_I.o setmcp_I.o setmix_I.o setsum_I.o setxuv_I.o setxv_I.o setxz_I.o solve_I.o \
     strsum_I.o xpot_I.o ypot_I.o \
     cofpotmpi_I.o getaldmpi_I.o getoldmpi_I.o getscdmpi_I.o \
     maneigmpi_I.o matrixmpi_I.o newcompi_I.o \
     improvmpi_I.o  \
     scfmpi_I.o  itjpo_I.o ispar_I.o\
     setcslmpi_I.o setdbgmpi_I.o setlagmpi_I.o lodcsh2GG_I.o lodcslmpiGG_I.o \
     systemmem_I.o systemfreemem_I.o setmcp_mem_I.o read1_mem_I.o read2_mem_I.o read3_mem_I.o \
     read02_mem_I.o read03_mem_I.o \
     \
     consis.o csfwgt.o dacon.o dampck.o dampor.o defcor.o \
     dsubrs.o eigen.o endsum.o engoutgg.o estim.o fco.o gco.o getaldwt.o \
     getoldwt.o in.o lagcon.o maxarr.o newe.o orbout.o \
     orthor.o orthy.o out.o outbnd.o prtrsl.o prwf.o setcof.o setdbg.o \
     setham.o setmcp.o setmix.o setsum.o setxuv.o setxv.o setxz.o solve.o \
     strsum.o xpot.o ypot.o \
     cofpotmpi.o getaldmpi.o getoldmpi.o getscdmpi.o hmoutmpi.o  \
     improvmpi.o maneigmpi.o matrixmpi.o newcompi.o  rscfmpivu.o scfmpi.o \
     setcslmpi.o setdbgmpi.o setlagmpi.o lodcsh2GG.o lodcslmpiGG.o itjpo.o ispar.o \
     systemmem.o systemfreemem.o setmcp_mem.o read1_mem.o read2_mem.o read3_mem.o \
     read02_mem.o read03_mem.o

$(EXE): $(OBJS)
	$(FC_MPI) -o $@ $? $(FC_MPILD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC_MPI) -c $(FC_MPIFLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
