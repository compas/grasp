.SUFFIXES: .f90 .mod

EXE = rci_mpi
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ../../lib
MODDIR = ${SRCLIBDIR}/libmod
MODL92 = ${SRCLIBDIR}/lib9290
MODLRANG90 = ${SRCLIBDIR}/librang90
MODLMCP90 = ${SRCLIBDIR}/libmcp90
MODLDVD = ${SRCLIBDIR}/libdvd90
MODLMPIU90 = ${SRCLIBDIR}/mpi90

GRASPLIBS =-l9290 -lmod  -lrang90 -lmcp90 -ldvd90 -lmpiu90 -l9290


APP_LIBS = -L ${GRASPLIB} ${GRASPLIBS} 

#   Define data types
VASTO = ${MODDIR}/vast_kind_param_M.o

APP_OBJ= \
        iabint_I.o	   funk_I.o     ncharg_I.o     skint_I.o \
	auxblk_I.o         funl_I.o     qed_I.o        snrc_I.o \
	bessel_I.o         fzalf_I.o    qed_slfen_I.o  spodmv_I.o \
	breid_I.o          genmat2_I.o  ratden_I.o     strsum_I.o \
	dspevx_I.o         genmat_I.o   rkint_I.o      talk_I.o \
        brint1_I.o         brint2_I.o   brint3_I.o     brint4_I.o \
        brint5_I.o  	   brint6_I.o   keint_I.o      vint_I.o \
        genintrkwrap_I.o   genintrk_I.o \
	vpint_I.o 	   lodmixmpi_I.o	getcid_I.o\
	brintf_I.o         hmout_I.o    rkintc_I.o     triangrk_I.o \
	brra_I.o           hovlap_I.o   vac2_I.o \
	cxk_I.o            indtpi_I.o   setdbg_I.o     vac4_I.o \
	klamaq_I.o   setmixmpi_I.o     vacpol_I.o \
	dnicmv_I.o         lodres_I.o   setres_I.o     vinti_I.o \
	maneig_I.o   setsum_I.o     zkf_I.o \
	engout_I.o         mohr_I.o     setham_gg_I.o\
	vpintf_I.o	   iniestdm_I.o iniestsd_I.o	wghtd5_I.o \
        matrix_I.o genintbreit1_I.o triangbreit1_I.o genintbreit2_I.o \
        triangbreit2_I.o genintbreit1wrap_I.o genintbreit2wrap_I.o \
        auxblk.o bessel.o breid.o brint1.o brint2.o brint3.o brint4.o \
        brint5.o brint6.o brintf.o brra.o cxk.o engout.o evcout.o funk.o \
        funl.o fzalf.o genintrk.o genmat.o hmout.o hovlap.o iabint.o indtpi.o \
        keint.o klamaq.o mohr.o ncharg.o qed.o rkint.o rkintc.o setdbg.o \
        setsum.o skint.o snrc.o talk.o triangrk.o vac2.o vac4.o vacpol.o \
        genintrkwrap.o vint.o vinti.o vpint.o vpintf.o wghtd5.o zkf.o \
        dnicmv.o genmat2.o getcid.o iniestdm.o iniestsd.o lodmixmpi.o lodres.o \
        maneig.o  rci90mpi.o  setham_gg.o setmixmpi.o setres.o \
        spodmv.o strsum.o qed_slfen.o matrix.o genintbreit1.o triangbreit1.o \
        genintbreit2.o triangbreit2.o genintbreit1wrap.o genintbreit2wrap.o

$(EXE): $(APP_OBJ)
	$(FC_MPI) -o $(BINFILE) $(FC_MPILD) $(APP_OBJ) $(APP_LIBS) \
		$(APP_LIBS)  $(LAPACK_LIBS) 

.f90.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) $< -I $(MODDIR) -I ${MODL92} -I $(MODLRANG90) -I $(MODLMCP90) \
	        -I $(MODDIR) -I $(MODLDVD) -I $(MODLMPIU90) -o $@

.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core  *.mod

