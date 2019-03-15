.SUFFIXES: .f90 .mod

LIBDIR = $(GRASP)/lib
LIBA = $(LIBDIR)/lib9290.a

LIBMOD = -lmod
LIBOBJ = $(OBJ)
MODDIR = ../libmod

RM = /bin/rm -f

MODOBJ = \


OBJ = \
   alcbuf_I.o calen_I.o clrx_I.o convrt_double_I.o convrt_I.o cord_I.o \
   cre_I.o cslh_I.o dinit_I.o dmerge_I.o factt_I.o convrt2_I.o \
   getyn_I.o ichkq1_I.o ichkq2_I.o ichop_I.o icopy_I.o iq_I.o  intrpq_I.o \
   ispar_I.o itjpo_I.o itrig_I.o jcup_I.o jqs_I.o   \
   lodcsh_I.o lodcsh2_I.o lodcsl_I.o lodiso_I.o ltab_I.o openfl_I.o orthsc_I.o\
   quad_I.o pack_I.o posfile_I.o interp_I.o\
   parsjl_I.o prsrcn_I.o prsrsl_I.o rint_I.o speak_I.o\
   setcsll_I.o setiso_I.o setj_I.o setmc_I.o setcsh_I.o setcsla_I.o setqna_I.o \
   skrc_I.o starttime_I.o stoptime_I.o count_I.o dpbdt_I.o draw_I.o rinti_I.o \
   slater_I.o yzk_I.o dcbsrw_I.o setpot_I.o cgamma_I.o arctan_I.o start_I.o \
   lodstate_I.o getrsl_I.o items_I.o ldigit_I.o setqic_I.o radgrd_I.o \
   nucpot_I.o setrwfa_I.o es_I.o lodrwf_I.o spicmv2_I.o iniest2_I.o \
   setcon_I.o \
   \
   alcbuf.o calen.o clrx.o convrt_double.o convrt.o cord.o \
   cre.o cslh.o dinit.o dmerge.o factt.o convrt2.o \
   getyn.o ichkq1.o ichkq2.o ichop.o icopy.o iq.o intrpq.o \
   ispar.o itjpo.o itrig.o jcup.o jqs.o  \
   lodcsh.o lodcsh2.o lodcsl.o  lodiso.o ltab.o openfl.o orthsc.o \
   quad.o pack.o posfile.o interp.o \
   parsjl.o prsrcn.o prsrsl.o rint.o speak.o\
   setcsll.o setiso.o setj.o setmc.o setcsh.o setcsla.o setqna.o \
   skrc.o starttime.o stoptime.o count.o dpbdt.o draw.o rinti.o \
   slater.o yzk.o dcbsrw.o setpot.o cgamma.o arctan.o start.o \
   lodstate.o getrsl.o items.o ldigit.o setqic.o radgrd.o \
   nucpot.o setrwfa.o es.o lodrwf.o spicmv2.o iniest2.o \
   setcon.o

install : $(LIBA)
$(LIBA) : $(LIBOBJ)
	@echo " Building " $(LIBA)
	ar -crvs $(LIBA)  $?

clean:
	-rm -f *.o *.mod core

.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I  $(MODDIR) -I .  -o $@
