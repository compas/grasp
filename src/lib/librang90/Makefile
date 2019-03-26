.SUFFIXES: .f90 .mod

LIBDIR = $(GRASP)/lib
LIBA = $(LIBDIR)/librang90.a

LIBMOD = -lmod
LIBOBJ = $(OBJ)
MODDIR = ../libmod
MODL92 = ../lib9290
MODLMCP90 = ../libmcp90

RM= /bin/rm -f

MODOBJ = \


OBJ = \
   a1jj_I.o awp1_I.o w1jjg_I.o awp1jjg_I.o wap1_I.o wap1jjg_I.o wj1_I.o \
   ww1_I.o suwjj_I.o nmtejj_I.o jfaze_I.o perko1_I.o perko2_I.o \
   itrexg_I.o itjj_I.o itjj2_I.o itjj3_I.o izas1_I.o mes_I.o eile_I.o \
   Rmeajj_I.o Rmeajj9_I.o  Rmeajj11_I.o  Rwjj_I.o \
   Rmew1jj_I.o Rmew3jj_I.o Rmew5jj_I.o Rmew7jj_I.o Rmew7bjj_I.o \
   c0t5s_I.o c1e0sm_I.o cle0sm_I.o c1e1sm_I.o jthn_I.o rumtjj_I.o Gracah1_I.o \
   dracah_I.o ittk_I.o nine0_I.o sixj1_I.o sixj35_I.o sixj4_I.o sixj_I.o \
   ixjtik_I.o nine_I.o sixj2_I.o sixj3_I.o sixj5_I.o snrc_I.o\
   onescalar_I.o onescalar1_I.o onescalar2_I.o \
   oneparticlejj1_I.o oneparticlejj2_I.o oneparticlejj_I.o \
   el1_I.o el2_I.o el3_I.o el31_I.o el32_I.o el33_I.o el4_I.o el41_I.o \
   el5_I.o el51_I.o el52_I.o el53_I.o rkco_gg_I.o \
   coulom_I.o diaga1_I.o diaga2_I.o diaga3_I.o diaga4_I.o diaga5_I.o \
   reco_I.o recop00_I.o recop1_I.o recop2_I.o rec3_I.o reco3_I.o \
   reco2_I.o reco4_I.o recoonescalar_I.o \
   gg12_I.o gg1112_I.o gg1122_I.o gg1222_I.o gg1233_I.o gg1234_I.o \
   \
   a1jj.o awp1.o w1jjg.o awp1jjg.o wap1.o wap1jjg.o wj1.o \
   ww1.o suwjj.o nmtejj.o jfaze.o perko1.o perko2.o \
   itrexg.o itjj.o itjj2.o itjj3.o izas1.o mes.o eile.o \
   Rmeajj.o Rmeajj9.o Rmeajj11.o Rwjj.o \
   Rmew1jj.o Rmew3jj.o Rmew5jj.o Rmew7jj.o Rmew7bjj.o \
   c0t5s.o c1e0sm.o cle0sm.o c1e1sm.o jthn.o rumtjj.o Gracah1.o \
   dracah.o ittk.o nine0.o sixj1.o sixj35.o sixj4.o sixj.o \
   ixjtik.o nine.o sixj2.o sixj3.o sixj5.o snrc.o \
   onescalar.o onescalar1.o onescalar2.o \
   oneparticlejj1.o oneparticlejj2.o oneparticlejj.o \
   el1.o el2.o el3.o el31.o el32.o el33.o el4.o el41.o \
   el5.o el51.o el52.o el53.o rkco_gg.o \
   coulom.o diaga1.o diaga2.o diaga3.o diaga4.o diaga5.o \
   reco.o recop00.o recop1.o recop2.o rec3.o reco3.o \
   reco2.o reco4.o recoonescalar.o \
   gg12.o gg1112.o gg1122.o gg1222.o gg1233.o gg1234.o

install : $(LIBA)
$(LIBA) : $(LIBOBJ)
	@echo " Building " $(LIBA)
	ar -crvs $(LIBA)  $?

.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I $(MODDIR) -I $(MODL92) -I $(MODLMCP90) -I $(MODDIR) -I .  -o $@

clean:
	-rm -f *.o *.mod core
