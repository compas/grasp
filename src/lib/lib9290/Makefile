LIBA=${GRASP}/lib/lib9290.a
MODULES_INSTALL=${GRASP}/lib/9290
FC_MODULES= -I ${GRASP}/src/lib/libmod

OBJS= \
	alcbuf.o \
	alcbuf_I.o \
	arctan.o \
	arctan_I.o \
	calen.o \
	calen_I.o \
	cgamma.o \
	cgamma_I.o \
	clrx.o \
	clrx_I.o \
	convrt.o \
	convrt_I.o \
	convrt2.o \
	convrt2_I.o \
	convrt_double.o \
	convrt_double_I.o \
	speak.o \
	speak_I.o \
	cord.o \
	cord_I.o \
	count.o \
	count_I.o \
	cre.o \
	cre_I.o \
	openfl.o \
	openfl_I.o \
	setcsll.o \
	setcsll_I.o \
	prsrsl.o \
	prsrsl_I.o \
	prsrcn.o \
	prsrcn_I.o \
	lodcsh.o \
	lodcsh_I.o \
	cslh.o \
	cslh_I.o \
	dcbsrw.o \
	dcbsrw_I.o \
	dinit.o \
	dinit_I.o \
	dmerge.o \
	dmerge_I.o \
	dpbdt.o \
	dpbdt_I.o \
	draw.o \
	draw_I.o \
	es.o \
	es_I.o \
	factt.o \
	factt_I.o \
	ldigit.o \
	ldigit_I.o \
	getrsl.o \
	getrsl_I.o \
	getyn.o \
	getyn_I.o \
	iq.o \
	iq_I.o \
	ichkq1.o \
	ichkq1_I.o \
	ichkq2.o \
	ichkq2_I.o \
	ichop.o \
	ichop_I.o \
	icopy.o \
	icopy_I.o \
	iniest2.o \
	iniest2_I.o \
	interp.o \
	interp_I.o \
	quad.o \
	quad_I.o \
	rint.o \
	rint_I.o \
	intrpq.o \
	intrpq_I.o \
	ispar.o \
	ispar_I.o \
	items.o \
	items_I.o \
	itjpo.o \
	itjpo_I.o \
	itrig.o \
	itrig_I.o \
	jcup.o \
	jcup_I.o \
	jqs.o \
	jqs_I.o \
	parsjl.o \
	parsjl_I.o \
	pack.o \
	pack_I.o \
	lodcsh2.o \
	lodcsh2_I.o \
	lodcsl.o \
	lodcsl_I.o \
	lodiso.o \
	lodiso_I.o \
	orthsc.o \
	orthsc_I.o \
	lodrwf.o \
	lodrwf_I.o \
	lodstate.o \
	lodstate_I.o \
	ltab.o \
	ltab_I.o \
	nucpot.o \
	nucpot_I.o \
	posfile.o \
	posfile_I.o \
	radgrd.o \
	radgrd_I.o \
	rinti.o \
	rinti_I.o \
	setcon.o \
	setcon_I.o \
	setcsh.o \
	setcsh_I.o \
	setcsla.o \
	setcsla_I.o \
	setiso.o \
	setiso_I.o \
	setj.o \
	setj_I.o \
	setmc.o \
	setmc_I.o \
	setpot.o \
	setpot_I.o \
	setqic.o \
	setqic_I.o \
	setqna.o \
	setqna_I.o \
	setrwfa.o \
	setrwfa_I.o \
	skrc.o \
	skrc_I.o \
	yzk.o \
	yzk_I.o \
	slater.o \
	slater_I.o \
	spicmv2.o \
	spicmv2_I.o \
	start.o \
	start_I.o \
	starttime.o \
	starttime_I.o \
	stoptime.o \
	stoptime_I.o

PHONY: install
install: $(LIBA)
	mkdir -p $(MODULES_INSTALL)
	cp -v *.mod $(MODULES_INSTALL)

$(LIBA): $(OBJS)
	@echo "Installing $@"
	ar -curs $@ $?

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(LIBA)
	-@rm *.o *.mod
	-@rm -R $(MODULES_INSTALL)
