LIBA=${GRASP}/lib/libmod.a

OBJS= \
	parameter_def_M.o \
	vast_kind_param_M.o \
	memory_man.o \
	AME_C.o \
	bcore_C.o \
	bess_C.o \
	bilst_C.o \
	biorb_C.o \
	blim_C.o \
	blk_C.o \
	blkidx_C.o \
	buffer_C.o \
	cffmat_C.o \
	ciimat_C.o \
	cimat_C.o \
	cnc_C.o \
	coeils_C.o \
	cons_C.o \
	core_C.o \
	corre_C.o \
	coun_C.o \
	couple_C.o \
	cteilsrk_C.o \
	cuto_C.o \
	damp_C.o \
	debug_C.o \
	decide_C.o \
	def_C.o \
	default_C.o \
	dumx_C.o \
	eigv_C.o \
	eigvec1_C.o \
	facts_C.o \
	fixd_C.o \
	foparm_C.o \
	fposition_C.o \
	grid_C.o \
	hblock_C.o \
	hmat_C.o \
	horb_C.o \
	hydpar_C.o \
	iccu_C.o \
	int_C.o \
	invt_C.o \
	iounit_C.o \
	jj2lsj_C.o \
	jj2lsjbio_C.o \
	jlabl_C.o \
	jqjc_C.o \
	keilst_C.o \
	kkstart_C.o \
	kkstartbreit_C.o \
	l1_C.o \
	l2_C.o \
	lagr_C.o \
	left_C.o \
	lib92p_C.o \
	lic13_C.o \
	m_C.o \
	mcp_C.o \
	mcpa_C.o \
	mcpb_C.o \
	mcpdata_C.o \
	mtjj2_C.o \
	mtjj_C.o \
	ncc_C.o \
	ncdist_C.o \
	node_C.o \
	npar_C.o \
	npot_C.o \
	offd_C.o \
	orb_C.o \
	orba_C.o \
	orbord_C.o \
	orthct_C.o \
	osc_C.o \
	ovl_C.o \
	peav_C.o \
	pos_C.o \
	pote_C.o \
	prnt_C.o \
	qedcut_C.o \
	rang_Int_C.o \
	ribojj11_C.o \
	ribojj9_C.o \
	ribojj_C.o \
	sacoef_C.o \
	sbc_C.o \
	sbdat1_C.o \
	sbdat_C.o \
	scf_C.o \
	stat_C.o \
	stor_C.o \
	syma_C.o \
	tatb_C.o \
	terms_C.o \
	titl_C.o \
	trk_C.o \
	vinlst_C.o \
	vpilst_C.o \
	wave_C.o \
	wchblk_C.o \
	wfac_C.o \
	where_C.o \
	whfrom_C.o

$(LIBA): $(OBJS)
	@echo "Installing $@"
	ar -curs $@ $?

%.o: %.f90
	$(FC) $(FC_FLAGS) -c -o $@ $<

clean:
	-rm -f $(LIBA)
	-rm -f *.o *.mod
