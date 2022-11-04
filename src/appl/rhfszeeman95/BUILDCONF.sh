EXE=rhfszeeman95
LIBRARIES="9290 mod rang90 mcp90 9290 mod"
LAPACK=true
FILES="
opt6_C.f90
engouth_I.f90
gethfd_I.f90
hfszeeman_I.f90
matelt_I.f90
rinthf_I.f90
setdbg_I.f90
setsum_I.f90
strsum_I.f90
getmixblock_I.f90
engouth.f90
gethfd.f90
hfszeeman.f90
hfszeeman06.f90
matelt.f90
rinthf.f90
setdbg.f90
setsum.f90
strsum.f90
getmixblock.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
