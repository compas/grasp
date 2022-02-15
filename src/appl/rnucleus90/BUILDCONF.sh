EXE=rnucleus
LIBRARIES="9290 mod"
FILES="
nucleus_m.f90
skfun.f90 skfun_I.f90
estrms.f90 estrms_I.f90
getcpr.f90 getcpr_I.f90
geniso.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
