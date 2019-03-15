.SUFFIXES: .f90 .mod

LIBDIR = $(GRASP)/lib
LIBA = $(LIBDIR)/libdvd90.a

LIBOBJ = $(OBJ)
MODDIR = ../libmod
MODL92 = ../lib9290
MODLMPIU90 = ../mpi90

RM= /bin/rm -f


OBJ = \
  adds_I.o       gdvd_I.o        mgs_nrm_I.o     ovflow_I.o     \
  dvdrvr_I.o     iniest_I.o      multbc_I.o      tstsel_I.o     \
  dvdson_I.o     initdvd_I.o     newvec_I.o     \
\
   dvdson.o      gdvd.o         iniest.o


install: $(LIBA)
$(LIBA) : $(LIBOBJ)
	@echo " Building " $(LIBA)
	ar -curvs $(LIBA) $(LIBOBJ)
	ranlib $(LIBA)

.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I $(MODDIR) -I $(MODL92) -I $(MODLMPIU90) -I .  -o $@

clean:
	rm -f *.o *.mod core
