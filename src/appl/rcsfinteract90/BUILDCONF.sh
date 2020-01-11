EXE=rcsfinteract
LIBRARIES="rang90 9290 mod"
FILES="
onescalar1INT.f90 onescalar1INT_I.f90
onescalar2INT.f90 onescalar2INT_I.f90
onescalarINT.f90 onescalarINT_I.f90
el1INT.f90 el1INT_I.f90
el2INT.f90 el2INT_I.f90
el31INT.f90 el31INT_I.f90
el32INT.f90 el32INT_I.f90
el33INT.f90 el33INT_I.f90
el3INT.f90 el3INT_I.f90
el41INT.f90 el41INT_I.f90
el4INT.f90 el4INT_I.f90
el51INT.f90 el51INT_I.f90
el52INT.f90 el52INT_I.f90
el53INT.f90 el53INT_I.f90
el5INT.f90 el5INT_I.f90
Interact_MR.f90 Interact_MR_I.f90
Interact_csf.f90 Interact_csf_I.f90
getinf.f90 getinf_I.f90
lodcsl_CSF.f90 lodcsl_CSF_I.f90
lodcsl_MR.f90 lodcsl_MR_I.f90
set_CSF_list.f90 set_CSF_list_I.f90
set_CSF_number.f90 set_CSF_number_I.f90

RCSFinteract.f90

# Not referenced in the original makefile
#recoonescalar.f90 recoonescalar_I.f90
#el52_I.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
