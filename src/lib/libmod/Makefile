.SUFFIXES: .f90 .mod

LIBDIR = $(GRASP)/lib
LIBA = $(LIBDIR)/libmod.a
F90 = $(FC)
F90_FLAGS = $(FC_FLAGS)

RM = /bin/rm -f

VASTO = vast_kind_param_M.o
PARDEF = parameter_def_M.o

MODOBJ = \
   memory_man.o \
   biorb_C.o  blk_C.o buffer_C.o cffmat_C.o  ciimat_C.o  cimat_C.o \
   cons_C.o couple_C.o cuto_C.o debug_C.o default_C.o hmat_C.o hblock_C.o iccu_C.o \
   \
   def_C.o default_C.o dumx_C.o grid_C.o iounit_C.o jj2lsj_C.o l1_C.o l2_C.o \
   \
   eigv_C.o facts_C.o foparm_C.o jqjc_C.o lib92p_C.o m_C.o mcp_C.o mtjj_C.o mtjj2_C.o \
   ncc_C.o npar_C.o orb_C.o orbord_C.o osc_C.o prnt_C.o left_C.o  eigvec1_C.o \
   jj2lsjbio_C.o mcpdata_C.o offd_C.o \
   ribojj_C.o ribojj9_C.o ribojj11_C.o sbc_C.o sbdat_C.o sbdat1_C.o stat_C.o syma_C.o \
\
   tatb_C.o terms_C.o titl_C.o trk_C.o hydpar_C.o qedcut_C.o\
\
   rang_Int_C.o AME_C.o wave_C.o int_C.o scf_C.o pos_C.o jlabl_C.o npot_C.o pote_C.o \
   damp_C.o invt_C.o coun_C.o lic13_C.o peav_C.o lagr_C.o fixd_C.o \
   ovl_C.o orba_C.o core_C.o corre_C.o cnc_C.o node_C.o wfac_C.o orthct_C.o sbc_C.o \
   where_C.o wchblk_C.o  whfrom_C.o  blkidx_C.o mcpa_C.o mcpb_C.o \
\
   decide_C.o ncdist_C.o coeils_C.o bilst_C.o keilst_C.o vinlst_C.o vpilst_C.o stor_C.o bess_C.o \
   bcore_C.o cteilsrk_C.o kkstart_C.o horb_C.o blim_C.o fposition_C.o sacoef_C.o \
   kkstartbreit_C.o  #   mpi_C.o   delete  for serial code only, delte mpi_C.i  


install : $(LIBA)
$(LIBA) : ${VASTO} ${PARDEF} $(MODOBJ)
	@echo " Building install libmod.a"
	ar -curs $(LIBA)  $?

${MODOBJ} : ${VASTO} ${PARDEF}

 clean:
	-rm -f *.o *.mod

.f90.o:
	@echo "Compiling source file $< into $@"
	$(F90) -c $(F90_FLAGS) $<

.f90.mod:
	$(F90) -c $(F90_FLAGS) $?

SRC =  \
   memory_man.f90 \
   blk_C.f90 buffer_C.f90 \
   cons_C.f90 couple_C.f90 debug_C.f90 default_C.f90 hmat_C.f90 hblock_C.f90 \
   hmat_C.f90 iccu_C.f90 \
   \
   def_C.f90 default_C.f90 dumx_C.f90 grid_C.f90 iounit_C.f90 jj2lsj_C.f90 l1_C.f90 l2_C.f90 \
   \
   eigv_C.f90 facts_C.f90 foparm_C.f90 lib92p_C.f90 m_C.f90 mcp_C.f90 mtjj_C.f90 mtjj2_C.f90 \
   mpi_C.f90 ncc_C.f90 npar_C.f90 orb_C.f90 prnt_C.f90 \
   \
   ribojj_C.f90 ribojj9_C.f90 ribojj11_C.f90 stat_C.f90 syma_C.f90 \
   itatb_C.f90 terms_C.f90 AME_C.f90 trk_C.f90 \
\
   rang_Int_C.f90 wave_C.f90 int_C.f90 scf_C.f90 pos_C.f90 jlabl_C.f90 npot_C.f90 pote_C.f90 \
   damp_C.f90 invt_C.f90 coun_C.f90 lic13_C.f90 peav_C.f90 lagr_C.f90 fixd_C.f90 \
   ovl_C.f90 orba_C.f90 core_C.f90 corre_C.f90 cnc_C.f90 node_C.f90 wfac_C.f90 orthct_C.f90 \
   sbc_C.f90 where_C.f90 wchblk_C.f90 blkidx_C.f90 mcpa_C.f90 mcpb_C.f90 \
\
   decide_C.f90 ncdist_C.f90 coeils_C.f90 bilst_C.f90 keilst_C.f90 vinlst_C.f90 vpilst_C.f90 stor_C.f90 bess_C.f90 \
   bcore_C.f90 cteilsrk_C.f90 kkstart_C.f90 horb_C.f90 qedcut_C.f90
