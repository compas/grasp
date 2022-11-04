EXE=rdensity
LIBRARIES="9290 mod rang90 mcp90 mod 9290" # Comment [JG, 2022]: Yes, this is strange - but mod and 9290 has to be linked twice like this in order for things to compile
LAPACK=true
FILES="
rdensity_C.f90
teilst_C.f90
cxk_I.f90
getmixblock_I.f90
getsmd_I.f90
polint_I.f90
rdensity_cal_I.f90
rintdens_I.f90
rintdensvec_I.f90
natorbnew_I.f90
setdbg_I.f90
setdens_I.f90
rdensity.f90
cxk.f90
getmixblock.f90
getsmd.f90
polint.f90
rdensity_cal.f90
rintdens.f90
rintdensvec.f90
natorbnew.f90
setdbg.f90
setdens.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
