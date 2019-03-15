.SUFFIXES: .f90 .mod

# executable :: jj2lsj90
EXE = jj2lsj
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ../../lib
MODDIR = ${SRCLIBDIR}/libmod
MODL92 = ${SRCLIBDIR}/lib9290
MODLRANG = ${SRCLIBDIR}/librang90
GRASPLIBS =  -l9290 -lrang90 -lmod

APP_LIBS = -L ${GRASPLIB} ${GRASPLIBS}

#   Define data types
VASTO = ${MODDIR}/vast_kind_param_M.o

#   Define Commons
Commons = jj2lsj_data_1_C.o jj2lsj_data_2_C.o jj2lsj_data_3_C.o

#  Define memory management module
Memory = ${MODDIR}/memory_man.o

#  Define interface to routines from the library
Interface = packLS_I.o getmixblock_I.o idigit_I.o lval_I.o 

APP_OBJ = \
   packLS.o getmixblock.o idigit.o lval.o \
   jj2lsj_code.o jj2lsj2K.o 

$(EXE): ${VASTO} ${Commons} ${Memory} ${Interface} $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(Commons) ${Interface} $(APP_OBJ) $(APP_LIBS)

.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I $(MODDIR) -I . -I ${MODL92} \
	-I ${MODLRANG} -o $@

.f90.mod:
	$(FC) -c $(FC_FLAGS) $< -I $(MODDIR) -I . -I ${MODL92} \
	-I ${MODLRANG} -o $@

clean:
	-rm -f *.o core *.mod

APP_SRC =  \
      jj2lsj_data_1.f90 jj2lsj_data_2.f90 jj2lsj_data_3.f90 \
      jj2lsj_code.f90 jj2lsj2K.f90
