# Make sure the test executables end up on build/test/, not build/bin/
# But we still need to know where the GRASP executables live, so we save that
# information first in a different variable.
set(BUILD_BINDIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY)

# Let's catch implicit uses of routines in the test prorams and libraries
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Werror=implicit-interface -Werror=implicit-procedure")

# Check lib92 routines -- this currently acts as a simple integration test,
# making sure that everything actually got compiled properly. The test checks
# that the QUAD routine from lib92 actually produces reasonable numbers.
add_executable(test.lib9290_quad
    lib9290_quad.f90
)
target_link_libraries_Fortran(test.lib9290_quad PUBLIC mod 9290)
add_test(lib9290_quad test.lib9290_quad)

# Simple integration test testing the rnucleus-rcsfgenerate-rwfnestimate-rangular-rmcdhf
# workflow.
foreach(conf "serial" "mpi" "mpi-longpath")
    add_test(NAME integration.workflow.${conf}
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/integration/mpitmp/run.sh "${conf}"
    )
    set_property(TEST integration.workflow.${conf} PROPERTY ENVIRONMENT
        "GRASP_BINDIR=${BUILD_BINDIR}"
    )
endforeach()
