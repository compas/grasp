install(PROGRAMS lscomp.pl DESTINATION bin/)
install(PROGRAMS rsave DESTINATION bin/)

add_executable(rasfsplit rasfsplit.f90)
target_link_libraries_Fortran(rasfsplit PRIVATE mod 9290)
install(TARGETS rasfsplit DESTINATION bin/)

add_executable(rcsfblock rcsfblock.f90)
target_link_libraries_Fortran(rcsfblock PRIVATE mod 9290)
install(TARGETS rcsfblock DESTINATION bin/)

add_executable(rcsfmr rcsfmr.f90)
target_link_libraries_Fortran(rcsfmr PRIVATE mod 9290)
install(TARGETS rcsfmr DESTINATION bin/)

add_executable(rcsfsplit rcsfsplit.f90)
target_link_libraries_Fortran(rcsfsplit PRIVATE mod 9290)
install(TARGETS rcsfsplit DESTINATION bin/)

add_executable(rhfs_lsj rhfs_lsj.f90)
target_link_libraries_Fortran(rhfs_lsj PRIVATE mod 9290)
install(TARGETS rhfs_lsj DESTINATION bin/)

add_executable(rlevelseV rlevelseV.f90)
target_link_libraries_Fortran(rlevelseV PRIVATE mod 9290)
install(TARGETS rlevelseV DESTINATION bin/)

add_executable(rlevels rlevels.f90)
target_link_libraries_Fortran(rlevels PRIVATE mod 9290)
install(TARGETS rlevels DESTINATION bin/)

add_executable(rmixaccumulate rmixaccumulate.f90)
target_link_libraries_Fortran(rmixaccumulate PRIVATE mod 9290)
install(TARGETS rmixaccumulate DESTINATION bin/)

add_executable(rmixextract rmixextract.f90)
target_link_libraries_Fortran(rmixextract PRIVATE mod 9290)
install(TARGETS rmixextract DESTINATION bin/)

add_executable(rseqenergy rseqenergy.f90)
target_link_libraries_Fortran(rseqenergy PRIVATE mod 9290)
install(TARGETS rseqenergy DESTINATION bin/)

add_executable(rseqhfs rseqhfs.f90)
target_link_libraries_Fortran(rseqhfs PRIVATE mod 9290)
install(TARGETS rseqhfs DESTINATION bin/)

add_executable(rseqtrans rseqtrans.f90)
target_link_libraries_Fortran(rseqtrans PRIVATE mod 9290)
install(TARGETS rseqtrans DESTINATION bin/)

add_executable(rtabhfs rtabhfs.f90)
target_link_libraries_Fortran(rtabhfs PRIVATE mod 9290)
install(TARGETS rtabhfs DESTINATION bin/)

add_executable(rtablevels rtablevels.f90)
target_link_libraries_Fortran(rtablevels PRIVATE mod 9290)
install(TARGETS rtablevels DESTINATION bin/)

add_executable(rtabtrans1 rtabtrans1.f90)
target_link_libraries_Fortran(rtabtrans1 PRIVATE mod 9290)
install(TARGETS rtabtrans1 DESTINATION bin/)

add_executable(rtabtrans2 rtabtrans2.f90)
target_link_libraries_Fortran(rtabtrans2 PRIVATE mod 9290)
install(TARGETS rtabtrans2 DESTINATION bin/)

add_executable(rtabtransE1 rtabtransE1.f90)
target_link_libraries_Fortran(rtabtransE1 PRIVATE mod 9290)
install(TARGETS rtabtransE1 DESTINATION bin/)

add_executable(rwfnmchfmcdf rwfnmchfmcdf.f90)
target_link_libraries_Fortran(rwfnmchfmcdf PRIVATE mod 9290)
install(TARGETS rwfnmchfmcdf DESTINATION bin/)

add_executable(rwfnplot rwfnplot.f90)
target_link_libraries_Fortran(rwfnplot PRIVATE mod 9290)
install(TARGETS rwfnplot DESTINATION bin/)

add_executable(rwfnrelabel rwfnrelabel.f90)
target_link_libraries_Fortran(rwfnrelabel PRIVATE mod 9290)
install(TARGETS rwfnrelabel DESTINATION bin/)

add_executable(rwfnrotate rwfnrotate.f90)
target_link_libraries_Fortran(rwfnrotate PRIVATE mod 9290)
install(TARGETS rwfnrotate DESTINATION bin/)

add_executable(wfnplot wfnplot.f90)
target_link_libraries_Fortran(wfnplot PRIVATE mod 9290)
install(TARGETS wfnplot DESTINATION bin/)

