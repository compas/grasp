EXE=jj2lsj
LIBRARIES="9290 rang90 mod "
FILES="
getmixblock.f90 getmixblock_I.f90
idigit.f90 idigit_I.f90
lval.f90 lval_I.f90
packLS.f90 packLS_I.f90
jj2lsj_data_1_C.f90 jj2lsj_data_2_C.f90 jj2lsj_data_3_C.f90
jj2lsj_code.f90
jj2lsj2K.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
