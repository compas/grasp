.SUFFIXES: .f90 .mod

EXE = rbiotransform
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ../../lib
MODDIR = ${SRCLIBDIR}/libmod
MODL92 = ${SRCLIBDIR}/lib9290
MODLRANG90 = ${SRCLIBDIR}/librang90
MODLMCP90 = ${SRCLIBDIR}/libmcp90
GRASPLIBS =-l9290 -lmod -lrang90 -lmcp90 -ldvd90

APP_LIBS = -L${GRASPLIB} ${GRASPLIBS}

APP_OBJ=  \
	angdata_I.o  getmix_I.o   invmat_I.o   pamtmt_I.o   setcslb_I.o \
	biotr1_I.o   gets_I.o     kapdata_I.o  prsym_I.o    setvec_I.o \
	bndinv_I.o   ichkq1_I.o   lodrwff_I.o  qqsort_I.o   tcsl_I.o \
	brkt_I.o     ielsum_I.o   lodrwfi_I.o  radfile_I.o  tiinig_I.o \
	citrag_I.o   ifnmnx_I.o   lulu_I.o     radpar_I.o   trpmat_I.o \
	copvec_I.o   inprod_I.o   matml4_I.o   rintff_I.o   ulla_I.o \
	fname_I.o    intrpqf_I.o  mcpout_gg_I.o   rintii_I.o   vecsum_I.o \
	genmcp_I.o   intrpqi_I.o  orbord_I.o   scalve_I.o   wrtmat_I.o \
        mcpin_I.o ti1tv_I.o lodcslBio_I.o \
	biotr.o fname.o intrpqf.o intrpqi.o kapdata.o lodrwff.o lodrwfi.o \
        radpar.o rintff.o rintii.o tcsl.o brkt.o gets.o wrtmat.o  \
        biotr1.o radfile.o ifnmnx.o ielsum.o copvec.o invmat.o \
        ulla.o trpmat.o matml4.o scalve.o setvec.o pamtmt.o bndinv.o lulu.o \
        inprod.o genmcp.o mcpin.o getmix.o qqsort.o citrag.o tiinig.o \
        ti1tv.o vecsum.o prsym.o setcslb.o angdata.o mcpout_gg.o lodcslBio.o

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) $(APP_LIBS) \
		$(APP_LIBS) ${LAPACK_LIBS} 

.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I $(MODDIR) -I . -I ${MODL92} -I ${MODLRANG90} -I $(MODLMCP90)  -o $@

.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core *.mod


