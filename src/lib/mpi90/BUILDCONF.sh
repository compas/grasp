LIB=mpiu90
LIBRARIES="mod 9290"
ISMPI=true
FILES="
mpi_C.f90 mpiu.f90

cpath.f90 cpath_I.f90
cslhmpi.f90 cslhmpi_I.f90
iniestmpi.f90 iniestmpi_I.f90
lodcslmpi.f90 lodcslmpi_I.f90
lodrwfmpi.f90 lodrwfmpi_I.f90
setisompi.f90 setisompi_I.f90
setrwfmpi.f90 setrwfmpi_I.f90
spicmvmpi.f90 spicmvmpi_I.f90
sys_chdir.f90 sys_chdir_I.f90
sys_getwd.f90 sys_getwd_I.f90
sys_mkdir.f90 sys_mkdir_I.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
