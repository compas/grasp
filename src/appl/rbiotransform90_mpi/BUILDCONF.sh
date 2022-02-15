EXE=rbiotransform_mpi
LIBRARIES="rang90 mpiu90 9290 mod"
LAPACK=true
ISMPI=true
FILES="
angdatampi.f90 angdatampi_I.f90
bndinv.f90 bndinv_I.f90
brkt.f90 brkt_I.f90
vecsum.f90 vecsum_I.f90
wrtmat.f90 wrtmat_I.f90
copvec.f90 copvec_I.f90
fname.f90 fname_I.f90
getmixmpi.f90 getmixmpi_I.f90
getsmpi.f90 getsmpi_I.f90
ielsum.f90 ielsum_I.f90
ifnmnx.f90 ifnmnx_I.f90
inprod.f90 inprod_I.f90
invmat.f90 invmat_I.f90
rintff.f90 rintff_I.f90
rintii.f90 rintii_I.f90
intrpqf.f90 intrpqf_I.f90
intrpqi.f90 intrpqi_I.f90
kapdata.f90 kapdata_I.f90
lodcslBio.f90 lodcslBio_I.f90
lodrwffmpi.f90 lodrwffmpi_I.f90
lodrwfimpi.f90 lodrwfimpi_I.f90
prsym.f90 prsym_I.f90
lulu.f90 lulu_I.f90
setvec.f90 setvec_I.f90
matml4.f90 matml4_I.f90
qqsortmpi.f90 qqsortmpi_I.f90
mcpinmpi.f90 mcpinmpi_I.f90
mcpoutmpi_gg.f90 mcpoutmpi_gg_I.f90
pamtmt.f90 pamtmt_I.f90
radfilempi.f90 radfilempi_I.f90
radparmpi.f90 radparmpi_I.f90
scalve.f90 scalve_I.f90
setcslampi.f90 setcslampi_I.f90
setcslbmpi.f90 setcslbmpi_I.f90
tcsl.f90 tcsl_I.f90
ti1tv.f90 ti1tv_I.f90
tiinigmpi.f90 tiinigmpi_I.f90
trpmat.f90 trpmat_I.f90
ulla.f90 ulla_I.f90
citragmpi.f90 citragmpi_I.f90
biotr1.f90 biotr1_I.f90
biotrmpi.f90

# These routines were not being compiled in the original makefile. The
# interface file for ORBORD this was still linked though.
#orbord.f90 orbord_I.f90
#tiinig.f90 tiinig_I.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
