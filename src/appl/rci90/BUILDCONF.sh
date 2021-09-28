EXE=rci
LIBRARIES="rang90 9290 dvd90 mod"
LAPACK=true
FILES="
# Missing interface
evcout.f90

# For these routines, only the interface files were linked in the original
# makefile.
#dspevx.f90 dspevx_I.f90
#ratden.f90 ratden_I.f90
dspevx_I.f90
ratden_I.f90

# Missing implementation, but linked in the original makefile
#dmerge_dnicmv_I.f90

# Missing implementation and not referenced in the original makefile
#inter_I.f90

# Not referenced in the original makefile
#setham_to_genmat2_C.f90
#t.f90

# Common block remnants
where_C.f90
ncdist_C.f90

talk.f90 talk_I.f90
cxk.f90 cxk_I.f90
bessel.f90 bessel_I.f90
dnicmv.f90 dnicmv_I.f90
engout.f90 engout_I.f90
funk.f90 funk_I.f90
funl.f90 funl_I.f90
mohr.f90 mohr_I.f90
klamaq.f90 klamaq_I.f90
fzalf.f90 fzalf_I.f90
breid.f90 breid_I.f90
zkf.f90 zkf_I.f90
rkint.f90 rkint_I.f90
rkintc.f90 rkintc_I.f90
skint.f90 skint_I.f90
brra.f90 brra_I.f90
brintf.f90 brintf_I.f90
brint1.f90 brint1_I.f90
brint2.f90 brint2_I.f90
brint3.f90 brint3_I.f90
brint4.f90 brint4_I.f90
brint5.f90 brint5_I.f90
brint6.f90 brint6_I.f90
triangbreit1.f90 triangbreit1_I.f90
triangbreit2.f90 triangbreit2_I.f90
genintbreit1.f90 genintbreit1_I.f90
genintbreit2.f90 genintbreit2_I.f90
genintrk.f90 genintrk_I.f90
getcid.f90 getcid_I.f90
hmout.f90 hmout_I.f90
hovlap.f90 hovlap_I.f90
iabint.f90 iabint_I.f90
indtpi.f90 indtpi_I.f90
iniestdm.f90 iniestdm_I.f90
iniestsd.f90 iniestsd_I.f90
keint.f90 keint_I.f90
lodmix.f90 lodmix_I.f90
lodres.f90 lodres_I.f90
spodmv.f90 spodmv_I.f90
maneig.f90 maneig_I.f90
ncharg.f90 ncharg_I.f90
qed.f90 qed_I.f90
qed_slfen.f90 qed_slfen_I.f90
setcsl.f90 setcsl_I.f90
setdbg.f90 setdbg_I.f90
setmix.f90 setmix_I.f90
setres.f90 setres_I.f90
setsum.f90 setsum_I.f90
shield.f90 shield_I.f90
strsum.f90 strsum_I.f90
triangrk.f90 triangrk_I.f90
vac2.f90 vac2_I.f90
vac4.f90 vac4_I.f90
vacpol.f90 vacpol_I.f90
vinti.f90 vinti_I.f90
vint.f90 vint_I.f90
vpintf.f90 vpintf_I.f90
vpint.f90 vpint_I.f90
wghtd5.f90 wghtd5_I.f90
setham_gg.f90 setham_gg_I.f90
genmat.f90 genmat_I.f90
genmat2.f90 genmat2_I.f90
auxblk.f90 auxblk_I.f90
matrix.f90 matrix_I.f90

rci92.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
