.SUFFIXES: .f90 .mod

LIBDIR = $(GRASP)/lib
LIBA = $(LIBDIR)/libmcp90.a

LIBMOD = -lmod
LIBOBJ = $(OBJ)
MODDIR = ../libmod
MODL92 = ../lib9290

RM = /bin/rm -f

MODOBJ = \


OBJ = \
   cxk_I.o talk_I.o \
   \
   cxk.o   talk.o

install : $(LIBA)
$(LIBA) : $(LIBOBJ)
	@echo " Building " $(LIBA)
	ar -crvs $(LIBA)  $?

clean:
	-rm -f *.o *.mod core

.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I  $(MODDIR) -I $(MODL92) -I .  -o $@
