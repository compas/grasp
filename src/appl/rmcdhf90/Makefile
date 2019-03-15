.SUFFIXES: .f90 .mod

EXE = rmcdhf
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ../../lib
MODDIR = ${SRCLIBDIR}/libmod
MODL92 = ${SRCLIBDIR}/lib9290
MODDVD = ${SRCLIBDIR}/libdvd90
GRASPLIBS =-l9290 -lmod -ldvd90

APP_LIBS = -L${GRASPLIB} ${GRASPLIBS}


APP_OBJ= mpi_s.o \
    itjpo_I.o ispar_I.o \
    dsubrs_I.o  getaldwt_I.o  orbout_I.o  setcof_I.o  setxuv_I.o\
    cofpot_I.o  eigen_I.o   getold_I.o    orthor_I.o  setcsl_I.o  setxv_I.o\
    consis_I.o  endsum_I.o  getoldwt_I.o  lagcon_I.o  orthy_I.o   setdbg_I.o  setxz_I.o\
    csfwgt_I.o  engoutgg_I.o  getscd_I.o    maneig_I.o  out_I.o     setham_I.o  solve_I.o\
    dacon_I.o   estim_I.o   hmout_I.o     matrix_I.o  outbnd_I.o  setlag_I.o  strsum_I.o\
    dampck_I.o  fco_I.o     improv_I.o    maxarr_I.o  prtrsl_I.o  setmcp_I.o  xpot_I.o\
    dampor_I.o  gco_I.o     in_I.o        newco_I.o   prwf_I.o    setmix_I.o  ypot_I.o\
    defcor_I.o  getald_I.o  newe_I.o      scf_I.o     setsum_I.o  lodcsh2GG_I.o \
    \
    itjpo.o ispar.o \
    rscfvu.o consis.o csfwgt.o dacon.o dampck.o dampor.o defcor.o \
    dsubrs.o eigen.o endsum.o engoutgg.o estim.o fco.o gco.o getaldwt.o \
    getoldwt.o in.o lagcon.o maxarr.o newe.o orbout.o \
    orthor.o orthy.o out.o outbnd.o prtrsl.o prwf.o setcof.o setdbg.o \
    setham.o setmcp.o setmix.o setsum.o setxuv.o setxv.o setxz.o solve.o \
    strsum.o xpot.o ypot.o \
    cofpot.o getald.o getold.o getscd.o hmout.o improv.o maneig.o \
    matrix.o newco.o scf.o setcsl.o setlag.o lodcsh2GG.o

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD)  $(APP_OBJ) $(APP_LIBS)  \
	${LAPACK_LIBS}

.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I $(MODDIR) -I ${MODL92} -I ${MODDVD} \
	        -o $@

.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o *.mod core
