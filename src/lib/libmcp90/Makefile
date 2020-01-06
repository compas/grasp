LIBA=${GRASP}/lib/libmcp.a
FC_MODULES=-I ${GRASP}/src/lib/libmod -I ${GRASP}/src/lib/lib9290 

OBJS= \
	cxk.o \
	cxk_I.o \
	talk.o \
	talk_I.o

$(LIBA): $(OBJS)
	@echo "Installing $@"
	ar -curs $@ $?

%.o: %.f90
	$(FC) $(FC_FLAGS) $(FC_MODULES) -c -o $@ $<

clean:
	-rm -f $(LIBA)
	-rm -f *.o *.mod
