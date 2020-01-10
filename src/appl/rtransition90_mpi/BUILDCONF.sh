EXE=rtransition_mpi
LIBRARIES="mcp90 dvd90 rang90 mpiu90 9290 mod"
LAPACK=true
ISMPI=true
FILES="
alclla.f90 alclla_I.f90
alcnma.f90 alcnma_I.f90
alcnsa.f90 alcnsa_I.f90
alcnta.f90 alcnta_I.f90
angdatampi.f90 angdatampi_I.f90
bessj.f90 bessj_I.f90
brkt.f90 brkt_I.f90
connect.f90 connect_I.f90
cpmix.f90 cpmix_I.f90
spme.f90 spme_I.f90
csfm.f90 csfm_I.f90
engout1.f90 engout1_I.f90
fname.f90 fname_I.f90
getrmpmpi.f90 getrmpmpi_I.f90
lodrwffmpi.f90 lodrwffmpi_I.f90
lodrwfimpi.f90 lodrwfimpi_I.f90
getosdmpi.f90 getosdmpi_I.f90
iqr.f90 iqr_I.f90
isparr.f90 isparr_I.f90
itjpor.f90 itjpor_I.f90
jcupr.f90 jcupr_I.f90
jqsr.f90 jqsr_I.f90
ldcsl1mpi.f90 ldcsl1mpi_I.f90
ldcsl2mpi.f90 ldcsl2mpi_I.f90
ldlbl1.f90 ldlbl1_I.f90
ldlbl2.f90 ldlbl2_I.f90
lodcslm.f90 lodcslm_I.f90
trsortmpi.f90 trsortmpi_I.f90
mctinmpi.f90 mctinmpi_I.f90
mctoutmpi_gg.f90 mctoutmpi_gg_I.f90
merg12mpi.f90 merg12mpi_I.f90
mrgcslmpi.f90 mrgcslmpi_I.f90
readmixmpi.f90 readmixmpi_I.f90
printaLS.f90 printaLS_I.f90
printa.f90 printa_I.f90
osclmpi.f90 osclmpi_I.f90
setcslm.f90 setcslm_I.f90
strsum.f90 strsum_I.f90
testmix.f90 testmix_I.f90
biosclmpi.f90

# Not referenced in the original makefile
#setcsl.f90 setcsl_I.f90
"
generate-makefile > Makefile
generate-cmakelists > CMakeLists.txt
