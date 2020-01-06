BIN=${GRASP}/bin
EXE=jj2lsj
LIBS=-lmod -l9290 -lrang90 

OBJS= \
	getmixblock.o \
	getmixblock_I.o \
	idigit.o \
	idigit_I.o \
	lval.o \
	lval_I.o \
	packLS.o \
	packLS_I.o \
	jj2lsj2K.o \
	jj2lsj_code.o \
	jj2lsj_data_1_C.o \
	jj2lsj_data_2_C.o \
	jj2lsj_data_3_C.o

$(BIN)/$(EXE): $(OBJS)
	$(FC) $(FC_LD) -o $@ $(LIBS) $?

%.o: %.f90
	$(FC) $(FC_FLAGS) -c -o $@ $<

clean:
	-rm -f *.o *.mod
