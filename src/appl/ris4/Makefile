EXE=${GRASP}/bin/ris4
LIBS=-L ${GRASP}/lib/ -lmod -l9290 -lrang90 -lmcp90 -lmod -l9290
FC_MODULES= -I ${GRASP}/src/lib/libmod -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/librang90 -I ${GRASP}/src/lib/libmcp90 -I ${GRASP}/src/lib/libmod -I ${GRASP}/src/lib/lib9290

OBJS= \
	ris_C.o \
	teilst_C.o \
	angdata_I.o \
	densread_I.o \
	densread_seltz_I.o \
	densnew_I.o \
	densnew_seltz_I.o \
	edensityfit_I.o \
	getmixblock_I.o \
	getsmd_I.o \
	polint_I.o \
	ris_cal_I.o \
	rintdens_I.o \
	rintdensvec_I.o \
	rinti_nms_I.o \
	rint_sms2_I.o \
	rint_sms3_I.o \
	setdbg_I.o \
	setsum_I.o \
	strsum_I.o \
	smsread_I.o \
	smsnew_I.o \
	sigma_1_I.o \
	sigma_2_I.o \
	vinti_I.o \
	engout_I.o \
	ris.o \
	angdata.o \
	densread.o \
	densread_seltz.o \
	densnew.o \
	densnew_seltz.o \
	edensityfit.o \
	getmixblock.o \
	getsmd.o \
	polint.o \
	ris_cal.o \
	rintdens.o \
	rintdensvec.o \
	rinti_nms.o \
	rint_sms2.o \
	rint_sms3.o \
	setdbg.o \
	setsum.o \
	strsum.o \
	smsread.o \
	smsnew.o \
	sigma_1.o \
	sigma_2.o \
	vinti.o \
	engout.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
