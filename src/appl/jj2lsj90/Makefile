BIN=${GRASP}/bin
EXE=jj2lsj
LIBS=-lmod -l9290 -lrang90 

OBJS= \
	getmixblock.f90 \
	getmixblock_I.f90 \
	idigit.f90 \
	idigit_I.f90 \
	lval.f90 \
	lval_I.f90 \
	packLS.f90 \
	packLS_I.f90 \
	jj2lsj2K.f90 \
	jj2lsj_code.f90 \
	jj2lsj_data_1_C.f90 \
	jj2lsj_data_2_C.f90 \
	jj2lsj_data_3_C.f90

$(BIN)/$(EXE): $(OBJS)
	$(FC) $(FC_LD) -o $@ $(LIBS) $?

%.o: %.f90
	$(FC) $(FC_FLAGS) -c -o $@ $<

clean:
	-rm -f *.o *.mod
