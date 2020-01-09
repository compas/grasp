EXE=rangular
LIBRARIES="rang90 mcp90 mpiu90 9290 mod"
LAPACK=true
ISMPI=true
FILES="
fndbeg.f90 fndbeg_I.f90
getinf.f90 getinf_I.f90
outsdampi.f90 outsdampi_I.f90
setdbg.f90 setdbg_I.f90
setdbgmpi.f90 setdbgmpi_I.f90
setmcp.f90 setmcp_I.f90
setmcpmpi.f90 setmcpmpi_I.f90
setsda.f90 setsda_I.f90
setsum.f90 setsum_I.f90
settmp.f90 settmp_I.f90
sort.f90 sort_I.f90
strsum.f90 strsum_I.f90
mcpmpi_gg.f90 mcpmpi_gg_I.f90
genmcpmpi.f90
"
generate-makefile > Makefile
generate-cmakelists > CMakeLists.txt
