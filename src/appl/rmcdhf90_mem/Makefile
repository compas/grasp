EXE=${GRASP}/bin/rmcdhf_mem
LIBS=-L ${GRASP}/lib/ -ldvd90 -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/libdvd90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

OBJS= mpi_s.o rmcdhf_mem_C.o \
    itjpo_I.o ispar_I.o \
    dsubrs_I.o  getaldwt_I.o  orbout_I.o  setcof_I.o  setxuv_I.o\
    cofpot_I.o  eigen_I.o   getold_I.o    orthor_I.o  setcsl_I.o  setxv_I.o\
    consis_I.o  endsum_I.o  getoldwt_I.o  lagcon_I.o  orthy_I.o   setdbg_I.o  setxz_I.o\
    csfwgt_I.o  engoutgg_I.o  getscd_I.o    maneig_I.o  out_I.o     setham_I.o  solve_I.o\
    dacon_I.o   estim_I.o   hmout_I.o     matrix_I.o  outbnd_I.o  setlag_I.o  strsum_I.o\
    dampck_I.o  fco_I.o     improv_I.o    maxarr_I.o  prtrsl_I.o  setmcp_I.o  xpot_I.o\
    dampor_I.o  gco_I.o     in_I.o        newco_I.o   prwf_I.o    setmix_I.o  ypot_I.o\
    defcor_I.o  getald_I.o  newe_I.o      scf_I.o     setsum_I.o  lodcsh2GG_I.o \
    systemmem_I.o systemfreemem_I.o setmcp_mem_I.o read1_mem_I.o read2_mem_I.o read3_mem_I.o \
    read02_mem_I.o read03_mem_I.o \
    \
    itjpo.o ispar.o \
    rscfvu.o consis.o csfwgt.o dacon.o dampck.o dampor.o defcor.o \
    dsubrs.o eigen.o endsum.o engoutgg.o estim.o fco.o gco.o getaldwt.o \
    getoldwt.o in.o lagcon.o maxarr.o newe.o orbout.o \
    orthor.o orthy.o out.o outbnd.o prtrsl.o prwf.o setcof.o setdbg.o \
    setham.o setmcp.o setmix.o setsum.o setxuv.o setxv.o setxz.o solve.o \
    strsum.o xpot.o ypot.o \
    cofpot.o getald.o getold.o getscd.o hmout.o improv.o maneig.o \
    matrix.o newco.o scf.o setcsl.o setlag.o lodcsh2GG.o \
    systemmem.o systemfreemem.o setmcp_mem.o read1_mem.o read2_mem.o read3_mem.o \
    read02_mem.o read03_mem.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
