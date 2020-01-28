EXE=${GRASP}/bin/hf

OBJS= \
	HF.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
