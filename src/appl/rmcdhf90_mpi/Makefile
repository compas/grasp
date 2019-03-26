.SUFFIXES: .f90 .mod

EXE = rmcdhf_mpi
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ../../lib
MODDIR = ${SRCLIBDIR}/libmod
MODL92 = ${SRCLIBDIR}/lib9290
MODDVD = ${SRCLIBDIR}/libdvd90
MODLMPIU90 = ${SRCLIBDIR}/mpi90
MODLMCP90 = ${SRCLIBDIR}/libmcp90
GRASPLIBS =-l9290 -lmod -lmcp90 -lmpiu90 -ldvd90
LAPACK90_LIB = -llapack -lblas

APP_LIBS = -L${GRASPLIB} ${GRASPLIBS}

#   Define data types
VASTO = ${MODDIR}/vast_kind_param_M.o

APP_OBJ= \
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
     \
     consis.o csfwgt.o dacon.o dampck.o dampor.o defcor.o \
     dsubrs.o eigen.o endsum.o engoutgg.o estim.o fco.o gco.o getaldwt.o \
     getoldwt.o in.o lagcon.o maxarr.o newe.o orbout.o \
     orthor.o orthy.o out.o outbnd.o prtrsl.o prwf.o setcof.o setdbg.o \
     setham.o setmcp.o setmix.o setsum.o setxuv.o setxv.o setxz.o solve.o \
     strsum.o xpot.o ypot.o \
     cofpotmpi.o getaldmpi.o getoldmpi.o getscdmpi.o hmoutmpi.o  \
     improvmpi.o maneigmpi.o matrixmpi.o newcompi.o  rscfmpivu.o scfmpi.o \
     setcslmpi.o setdbgmpi.o setlagmpi.o lodcsh2GG.o lodcslmpiGG.o itjpo.o ispar.o

$(EXE): $(APP_OBJ)
	$(FC_MPI) -o $(BINFILE) $(FC_MPILD) $(APP_OBJ) $(APP_LIBS) \
		$(APP_LIBS) ${LAPACK_LIBS}

.f90.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) $< -I $(MODDIR) -I ${MODL92} -I $(MODLMPIU90) -I ${MODDVD} \
                -o $@

.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o *.mod core
