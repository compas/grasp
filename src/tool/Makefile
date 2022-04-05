LIBS=-L ${GRASP}/lib/ -l9290 -lmod
FC_MODULES=-I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

all:  ${GRASP}/bin/lscomp.pl ${GRASP}/bin/rsave ${GRASP}/bin/rasfsplit ${GRASP}/bin/rcsfblock ${GRASP}/bin/rcsfmr ${GRASP}/bin/rcsfsplit ${GRASP}/bin/rhfs_lsj ${GRASP}/bin/rlevelseV ${GRASP}/bin/rlevels ${GRASP}/bin/rmixaccumulate ${GRASP}/bin/rmixextract ${GRASP}/bin/rseqenergy ${GRASP}/bin/rseqhfs ${GRASP}/bin/rseqtrans ${GRASP}/bin/rtabhfs ${GRASP}/bin/rtablevels ${GRASP}/bin/rtabtrans1 ${GRASP}/bin/rtabtrans2 ${GRASP}/bin/rtabtransE1 ${GRASP}/bin/rwfnmchfmcdf ${GRASP}/bin/rwfnplot ${GRASP}/bin/rwfnrelabel ${GRASP}/bin/rwfnrotate ${GRASP}/bin/wfnplot ${GRASP}/bin/rwfnpyplot ${GRASP}/bin/rwfntotxt

${GRASP}/bin/lscomp.pl: lscomp.pl
	cp $^ $@
	chmod u+x $@

${GRASP}/bin/rsave: rsave
	cp $^ $@
	chmod u+x $@

${GRASP}/bin/rasfsplit: rasfsplit.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rcsfblock: rcsfblock.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rcsfmr: rcsfmr.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rcsfsplit: rcsfsplit.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rhfs_lsj: rhfs_lsj.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rlevelseV: rlevelseV.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rlevels: rlevels.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rmixaccumulate: rmixaccumulate.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rmixextract: rmixextract.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rseqenergy: rseqenergy.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rseqhfs: rseqhfs.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rseqtrans: rseqtrans.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rtabhfs: rtabhfs.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rtablevels: rtablevels.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rtabtrans1: rtabtrans1.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rtabtrans2: rtabtrans2.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rtabtransE1: rtabtransE1.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rwfnmchfmcdf: rwfnmchfmcdf.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rwfnplot: rwfnplot.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rwfnrelabel: rwfnrelabel.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rwfnrotate: rwfnrotate.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/wfnplot: wfnplot.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rwfntotxt: rwfntotxt.o
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

${GRASP}/bin/rwfnpyplot: rwfnpyplot
	cp $^ $@
	chmod u+x $@

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-rm -f  ${GRASP}/bin/lscomp.pl ${GRASP}/bin/rsave ${GRASP}/bin/rasfsplit ${GRASP}/bin/rcsfblock ${GRASP}/bin/rcsfmr ${GRASP}/bin/rcsfsplit ${GRASP}/bin/rhfs_lsj ${GRASP}/bin/rlevelseV ${GRASP}/bin/rlevels ${GRASP}/bin/rmixaccumulate ${GRASP}/bin/rmixextract ${GRASP}/bin/rseqenergy ${GRASP}/bin/rseqhfs ${GRASP}/bin/rseqtrans ${GRASP}/bin/rtabhfs ${GRASP}/bin/rtablevels ${GRASP}/bin/rtabtrans1 ${GRASP}/bin/rtabtrans2 ${GRASP}/bin/rtabtransE1 ${GRASP}/bin/rwfnmchfmcdf ${GRASP}/bin/rwfnplot ${GRASP}/bin/rwfnrelabel ${GRASP}/bin/rwfnrotate ${GRASP}/bin/wfnplot ${GRASP}/bin/rwfntotxt ${GRASP}/bin/rwfnpyplot
	-rm -f *.o *.mod
