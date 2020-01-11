EXE=rangular
LIBRARIES="rang90 mcp90 9290 mod"
LAPACK=true
FILES="
allocCheck.f90 allocCheck_I.f90
fndbeg.f90 fndbeg_I.f90
getinf.f90 getinf_I.f90
sort.f90 sort_I.f90
sortmem.f90 sortmem_I.f90
outsda.f90 outsda_I.f90
setsda.f90 setsda_I.f90
mcp_gg.f90 mcp_gg_I.f90
setdbg.f90 setdbg_I.f90
setmcp.f90 setmcp_I.f90
setmcp2.f90 setmcp2_I.f90
setsum.f90 setsum_I.f90
settmpGG.f90 settmpGG_I.f90
strsum.f90 strsum_I.f90

genmcp.f90

# cons_C was not being compiled in the original makefile
#cons_C.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
