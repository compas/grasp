BIN=${GRASP}/bin
EXE=HF

OBJS= \
	HF.f90

$(BIN)/$(EXE): $(OBJS)
	$(FC) $(FC_LD) -o $@ $?

%.o: %.f90
	$(FC) $(FC_FLAGS) -c -o $@ $<

clean:
	-rm -f *.o *.mod
