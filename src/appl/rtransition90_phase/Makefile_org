.SUFFIXES: .f90 .mod

FC = gfortran                                                        # Fortran compiler 
FC_FLAGS = -O2 -fno-automatic                                        # Serial code compiler flags
LAPACK_LIBS = -llapack -lblas                                        # Lapack libraries

# Location of the 2018 root directory
GRASP = /Users/wenxianli/Documents/Coding/grasp2018_master_git_F95

EXE = rtransition_phase
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ${GRASP}/src/lib
MODDIR = ${SRCLIBDIR}/libmod
MODL92 = ${SRCLIBDIR}/lib9290
MODLRANG90 = ${SRCLIBDIR}/librang90
MODLMCP90 = ${SRCLIBDIR}/libmcp90
GRASPLIBS =-l9290 -lmod -lrang90 -lmcp90 -ldvd90

APP_LIBS = -L${GRASPLIB} ${GRASPLIBS}

APP_OBJ= \
	alclla_I.o   cpmix_I.o    itjpor_I.o   mctout_gg_I.o   spme_I.o\
	alcnma_I.o   engout1_I.o  jcupr_I.o    merg12_I.o   strsum_I.o\
	alcnsa_I.o   fname_I.o    jqsr_I.o     mrgcsl_I.o   testmix_I.o\
	alcnta_I.o   getosd_I.o   ldcsl2_I.o   oscl_I.o     trsort_I.o\
	angdata_I.o  getrmp_I.o   lodcslm_I.o  printa_I.o   ldcsl1_I.o \
	bessj_I.o    ichkq1_I.o   lodrwff_I.o  readmix_I.o  csfm_I.o\
	brkt_I.o     iqr_I.o      lodrwfi_I.o  setcsl_I.o\
	connect_I.o  isparr_I.o   mctin_I.o    setcslm_I.o\
	ldlbl1_I.o   ldlbl2_I.o   printaLS_I.o\
\
	bioscl.o fname.o ldcsl1.o ldcsl2.o lodcslm.o merg12.o mrgcsl.o \
         readmix.o iqr.o isparr.o itjpor.o jcupr.o jqsr.o lodrwfi.o \
         lodrwff.o getosd.o brkt.o getrmp.o strsum.o oscl.o connect.o \
         alcnsa.o alcnta.o  mctin.o csfm.o printa.o spme.o testmix.o trsort.o \
         bessj.o alclla.o alcnma.o setcslm.o angdata.o engout1.o mctout_gg.o \
         cpmix.o \
	ldlbl1.o ldlbl2.o printaLS.o

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) $(APP_LIBS) \
		$(APP_LIBS) ${LAPACK_LIBS} 

.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I $(MODDIR) -I . -I ${MODL92} -I ${MODLRANG90} -I $(MODLMCP90) -o $@

.f.o:
	$(F90) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core *.mod

