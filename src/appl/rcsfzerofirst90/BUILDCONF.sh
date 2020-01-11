EXE=rcsfzerofirst
LIBRARIES="9290 mod"
FILES="
lodcsl_Part.f90 lodcsl_Part_I.f90
lodcsl_Zero.f90 lodcsl_Zero_I.f90
set_CSF_number.f90 set_CSF_number_I.f90
set_CSF_ZFlist.f90 set_CSF_ZFlist_I.f90
RCSFzerofirst.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
