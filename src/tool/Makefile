.SUFFIXES: .f90 .mod

BIN = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
SRCLIBDIR = ../lib
MODDIR = ${SRCLIBDIR}/libmod
MODL9290 = ${SRCLIBDIR}/lib9290
GRASPLIBS = -l9290 -lmod 

APP_LIBS = -L${GRASPLIB} ${GRASPLIBS} -llapack -lblas

UTIL = rcsfsplit rmixaccumulate rseqenergy \
	rseqhfs rseqtrans rtablevels rtabtransE1 \

install: EXE 
	cp rsave $(GRASP)/bin
	cp lscomp.pl $(GRASP)/bin

EXE :        $(BIN)/rcsfsplit \
		$(BIN)/rmixaccumulate\
		$(BIN)/rseqenergy\
		$(BIN)/rseqhfs\
		$(BIN)/rseqtrans\
		$(BIN)/rtablevels\
		$(BIN)/rtabtransE1\
                $(BIN)/rcsfblock  \
		$(BIN)/rwfnmchfmcdf  \
                $(BIN)/rlevels    \
                $(BIN)/rtabtrans1  \
                $(BIN)/rwfnplot \
		$(BIN)/rtabtrans2\
		$(BIN)/rlevelseV\
		$(BIN)/rtabhfs\
		$(BIN)/wfnplot\
                $(BIN)/rmixextract \
		$(BIN)/rasfsplit\
                $(BIN)/rwfnrotate \
                $(BIN)/rwfnrelabel \
	        $(BIN)/rhfs_lsj    \
	        $(BIN)/rcsfmr    \
   

OBJ : rasfsplit.o rcsfsplit.o rmixaccumulate.o rseqenergy.o \
      rseqhfs.o rseqtrans.o rtablevels.o rtabtransE1.o \
      rmixextract.o rcsfblock.o  rwfnmchfmcdf.o \
      rlevels.o rtabtrans1.o format_mix.o rlevelsj.o \
      rwfnrotate.o  \
      rwfnplot.o rtabtrans2.o rlevelseV.o rtabhfs.o wfnplot.o rwfnrelabel.o \    
	rhfs_lsj.o rcsfmr.o\

$(BIN)/rasfsplit: rasfsplit.o
	$(FC) -o $(BIN)/rasfsplit  rasfsplit.f90 

$(BIN)/rcsfsplit: rcsfsplit.o
	$(FC) -o $(BIN)/rcsfsplit rcsfsplit.o 

$(BIN)/rmixaccumulate: rmixaccumulate.o
	$(FC) -o $(BIN)/rmixaccumulate rmixaccumulate.o 

$(BIN)/rseqenergy: rseqenergy.o
	$(FC) -o $(BIN)/rseqenergy rseqenergy.o 

$(BIN)/rseqtrans: rseqtrans.o
	$(FC) -o $(BIN)/rseqtrans rseqtrans.o 

$(BIN)/rseqhfs: rseqhfs.o
	$(FC) -o $(BIN)/rseqhfs rseqhfs.o 

$(BIN)/rtablevels: rtablevels.o
	$(FC) -o $(BIN)/rtablevels rtablevels.o 

$(BIN)/rtabtransE1: rtabtransE1.o
	$(FC) -o $(BIN)/rtabtransE1 rtabtransE1.o 

$(BIN)/format_mix  : format_mix.o
	$(FC) -o $(BIN)/format_mix  $(FC_LD) -L$(LIBDIR) \
                 format_mix.o  

$(BIN)/rmixextract: rmixextract.o
	$(FC) -o $(BIN)/rmixextract $(FC_LD)   rmixextract.o  $(APP_LIBS)

$(BIN)/rcsfblock: rcsfblock.o
	$(FC) -o $(BIN)/rcsfblock $(FC_LD) rcsfblock.o $(APP_LIBS)

$(BIN)/rwfnmchfmcdf: rwfnmchfmcdf.o
	$(FC) -o $(BIN)/rwfnmchfmcdf $(FC_LD) rwfnmchfmcdf.o $(APP_LIBS) 

$(BIN)/rlevels: rlevels.o
	$(FC) -o $(BIN)/rlevels $(FC_LD)  rlevels.o  $(APP_LIBS)

$(BIN)/rtabtrans1: rtabtrans1.o
	$(FC) -o $(BIN)/rtabtrans1 $(FC_LD)  rtabtrans1.o $(APP_LIBS) 

$(BIN)/rlevelsj: rlevelsj.o
	$(FC) -o $(BIN)/rlevelsj $(FC_LD) rlevelsj.o $(APP_LIBS) 

$(BIN)/rwfnplot: rwfnplot.o
	$(FC) -o $(BIN)/rwfnplot rwfnplot.o 

$(BIN)/rtabtrans2: rtabtrans2.o
	$(FC) -o $(BIN)/rtabtrans2 rtabtrans2.o 

$(BIN)/rlevelseV: rlevelseV.o
	$(FC) -o $(BIN)/rlevelseV $(FC_LD) rlevelseV.o $(APP_LIBS) 

$(BIN)/rtabhfs: rtabhfs.o
	$(FC) -o $(BIN)/rtabhfs rtabhfs.o 

$(BIN)/wfnplot: wfnplot.o
	$(FC) -o $(BIN)/wfnplot wfnplot.o 

$(BIN)/rwfnrotate: rwfnrotate.o
	$(FC) -o $(BIN)/rwfnrotate $(FC_LD) rwfnrotate.o $(APP_LIBS) 

$(BIN)/rwfnrelabel: rwfnrelabel.o
	$(FC) -o $(BIN)/rwfnrelabel $(FC_LD) rwfnrelabel.o $(APP_LIBS) 

$(BIN)/rhfs_lsj: rhfs_lsj.o
	$(FC) -o $(BIN)/rhfs_lsj  $(FC_LD) rhfs_lsj.o $(APPS_LIBS)

$(BIN)/rcsfmr: rcsfmr.o
	$(FC) -o $(BIN)/rcsfmr  $(FC_LD) rcsfmr.o $(APPS_LIBS)


.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I $(MODDIR) -I ${MODL9290} -I $(MODDIR) -o $@

clean:
	-rm -f *.o core *.mod



