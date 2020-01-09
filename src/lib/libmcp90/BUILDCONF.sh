LIB=mcp90
LIBRARIES="mod 9290"
FILES="
cxk.f90 cxk_I.f90
talk.f90 talk_I.f90
"
generate-makefile > Makefile
generate-cmakelists > CMakeLists.txt
