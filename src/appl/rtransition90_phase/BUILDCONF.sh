EXE=rtransition_phase
LIBRARIES="dvd90 rang90 mcp90 9290 mod"
LAPACK=true
FILES="
alclla.f90 alclla_I.f90
alcnma.f90 alcnma_I.f90
alcnsa.f90 alcnsa_I.f90
alcnta.f90 alcnta_I.f90
angdata.f90 angdata_I.f90
bessj.f90 bessj_I.f90
brkt.f90 brkt_I.f90
connect.f90 connect_I.f90
cpmix.f90 cpmix_I.f90
spme.f90 spme_I.f90
csfm.f90 csfm_I.f90
engout1.f90 engout1_I.f90
fname.f90 fname_I.f90
lodrwff.f90 lodrwff_I.f90
lodrwfi.f90 lodrwfi_I.f90
getrmp.f90 getrmp_I.f90
getosd.f90 getosd_I.f90
iqr.f90 iqr_I.f90
isparr.f90 isparr_I.f90
itjpor.f90 itjpor_I.f90
jcupr.f90 jcupr_I.f90
jqsr.f90 jqsr_I.f90
ldcsl1.f90 ldcsl1_I.f90
ldcsl2.f90 ldcsl2_I.f90
ldlbl1.f90 ldlbl1_I.f90
ldlbl2.f90 ldlbl2_I.f90
lodcslm.f90 lodcslm_I.f90
mctin.f90 mctin_I.f90
trsort.f90 trsort_I.f90
mctout_gg.f90 mctout_gg_I.f90
merg12.f90 merg12_I.f90
mrgcsl.f90 mrgcsl_I.f90
readmix.f90 readmix_I.f90
printaLS.f90 printaLS_I.f90
printa.f90 printa_I.f90
oscl.f90 oscl_I.f90
setcslm.f90 setcslm_I.f90
strsum.f90 strsum_I.f90
testmix.f90 testmix_I.f90
bioscl.f90

# Implementation not compiled in the original makefile
#ichkq1.f90 ichkq1_I.f90
#setcsl.f90 setcsl_I.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
