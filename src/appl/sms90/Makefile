EXE=${GRASP}/bin/sms
LIBS=-L ${GRASP}/lib/ -lrang90 -lmcp90 -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/librang90 -I ${GRASP}/src/lib/libmcp90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

OBJS= \
	dvpot_C.o \
	sms1_C.o \
	teilst_C.o \
	densmcp.o \
	densmcp_I.o \
	densnew.o \
	densnew_I.o \
	engout.o \
	engout_I.o \
	gco.o \
	gco_I.o \
	getmixblock.o \
	getmixblock_I.o \
	getsmd.o \
	getsmd_I.o \
	polint.o \
	polint_I.o \
	rintdens.o \
	rintdens_I.o \
	rintiso.o \
	rintiso_I.o \
	setdbg.o \
	setdbg_I.o \
	setmcp.o \
	setmcp_I.o \
	setsum.o \
	setsum_I.o \
	wghtd5.o \
	wghtd5_I.o \
	strsum.o \
	strsum_I.o \
	vinti.o \
	vinti_I.o \
	smsmcp.o \
	smsmcp_I.o \
	smsnew.o \
	smsnew_I.o \
	sms.o \
	sms_I.o \
	sms92.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-rm -f $(EXE)
	-rm -f *.o *.mod
