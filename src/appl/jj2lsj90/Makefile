EXE=${GRASP}/bin/jj2lsj
LIBS=-L ${GRASP}/lib/ -l9290 -lrang90 -lmod
FC_MODULES= -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/librang90 -I ${GRASP}/src/lib/libmod

OBJS= \
	getmixblock.o \
	getmixblock_I.o \
	idigit.o \
	idigit_I.o \
	lval.o \
	lval_I.o \
	packLS.o \
	packLS_I.o \
	jj2lsj_data_1_C.o \
	jj2lsj_data_2_C.o \
	jj2lsj_data_3_C.o \
	jj2lsj_code.o \
	jj2lsj2K.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-rm -f $(EXE)
	-rm -f *.o *.mod
