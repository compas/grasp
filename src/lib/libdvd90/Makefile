LIBA=${GRASP}/lib/libdvd90.a
FC_MODULES= -I ${GRASP}/src/lib/libmod

OBJS= \
	adds_I.o \
	dvdrvr_I.o \
	dvdson_I.o \
	initdvd_I.o \
	mgs_nrm_I.o \
	multbc_I.o \
	newvec_I.o \
	ovflow_I.o \
	tstsel_I.o \
	dvdson.o \
	gdvd.o \
	gdvd_I.o \
	iniest.o \
	iniest_I.o

$(LIBA): $(OBJS)
	@echo "Installing $@"
	ar -curs $@ $?

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-rm -f $(LIBA)
	-rm -f *.o *.mod
