add_library(mcp STATIC
    cxk.f90
    cxk_I.f90
    talk.f90
    talk_I.f90
)
setup_fortran_modules(mcp)
target_link_libraries_Fortran(mcp PRIVATE mod 9290)
