LIB=dvd90
LIBRARIES="mod"
FILES="
# dvdson.f90 contains the implementation for a bunch of routines, but we
# have separate interface files for them
adds_I.f90 dvdrvr_I.f90 dvdson_I.f90 initdvd_I.f90 mgs_nrm_I.f90
multbc_I.f90 newvec_I.f90 ovflow_I.f90 tstsel_I.f90
dvdson.f90

gdvd.f90 gdvd_I.f90
iniest.f90 iniest_I.f90
"
generate-makefile > Makefile
generate-cmakelists > CMakeLists.txt
