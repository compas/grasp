EXE=rwfnestimate
LIBRARIES="rang90 9290 mod"
LAPACK=true
FILES="
tail.f90 tail_I.f90
sbstep.f90 sbstep_I.f90
solvh.f90 solvh_I.f90
frmhyd.f90 frmhyd_I.f90
frmrwf.f90 frmrwf_I.f90
frmtfp.f90 frmtfp_I.f90
tfpot.f90 tfpot_I.f90
prtrem.f90 prtrem_I.f90
screenpar.f90 screenpar_I.f90
setdbg.f90 setdbg_I.f90
setsum.f90 setsum_I.f90
strsum.f90 strsum_I.f90
summry.f90 summry_I.f90
wrtrwf.f90 wrtrwf_I.f90
genrwf.f90 genrwf_I.f90
getinfo.f90 getinf_I.f90 # subroutine GETINF, implementation in getinfo.f90
erwf.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
