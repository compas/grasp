EXE=jjgen
LIBRARIES="mod 9290 rang90"
FILES="
lockad.f90 lockad_I.f90
reada.f90 reada_I.f90
lasa1.f90 lasa1_I.f90
lasa2.f90 lasa2_I.f90
adder.f90 adder_I.f90
slug.f90 slug_I.f90
kopp1.f90 kopp1_I.f90
kopp2.f90 kopp2_I.f90

# gen.f90 was not compiled in the Makefile, but genb.f90 was. genb.f90
# contains the implementation for GEN.
genb.f90 gen_I.f90

sluggo.f90 sluggo_I.f90
test.f90 test_I.f90
mergeb.f90 mergeb_I.f90
blanda.f90 blanda_I.f90
blandb.f90 blandb_I.f90
blandc.f90 blandc_I.f90
copy7t9.f90 copy7t9_I.f90
fivefirst.f90 fivefirst_I.f90
fivelines.f90 fivelines_I.f90
lika.f90 lika_I.f90
matain.f90 matain_I.f90
matbin.f90 matbin_I.f90
matcin.f90 matcin_I.f90
merge.f90 merge_I.f90
open79.f90 open79_I.f90
reffa.f90 reffa_I.f90

# There is also jjgen15.f90, but that was not being compiled in the original
# makefile.
jjgen15b.f90

# There are two more files that were not compiled in the original makefile
#lasax-reada.f90
#m.f90
"
generate-makefile > Makefile
generate-cmakelists > CMakeLists.txt
