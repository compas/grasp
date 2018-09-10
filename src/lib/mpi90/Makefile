.SUFFIXES: .f90 .mod

LIBDIR = $(GRASP)/lib
LIBA = $(LIBDIR)/libmpiu90.a

LIBOBJ = $(OBJ)
MODDIR = ../libmod
MODL92 = ../lib9290
PARDEF = $(MODDIR)/parameter_def_M.o

RM = /bin/rm -f

MODOBJ = \


OBJ =  mpi_C.o  \
   cpath_I.o cslhmpi_I.o iniestmpi_I.o lodcslmpi_I.o \
   lodrwfmpi_I.o setisompi_I.o setrwfmpi_I.o \
   spicmvmpi_I.o sys_chdir_I.o sys_getwd_I.o\
   sys_mkdir_I.o \
       \
   cpath.o cslhmpi.o iniestmpi.o lodcslmpi.o \
   lodrwfmpi.o      mpiu.o           setisompi.o      setrwfmpi.o \
   spicmvmpi.o      sys_chdir.o  sys_getwd.o  \
   sys_mkdir.o

install : $(LIBA)
$(LIBA) : $(LIBOBJ) ${PARDEF} $(MODOBJ)
	@echo " Building " $(LIBA)
	ar -crvs $(LIBA)  $?

clean :
	-rm -f *.o *.mod core

.f90.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) $<  -I $(MODDIR) -I $(MODL92) -I . -o $@
