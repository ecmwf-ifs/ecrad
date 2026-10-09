# (C) Copyright 2026- ECMWF.
# Licensed under the Apache Licence Version 2.0; see COPYING.

ecbuild_add_test(
    TARGET ecrad_${PREC}_flux_helper_ranges
    SOURCES flux_helper_ranges.F90
    LIBS ecrad.${PREC} ecrad_base.${PREC} OpenMP::OpenMP_Fortran
)

# NVHPC searches the module output directory before any -I path, so this must
# not inherit the last precision left by the top-level loop.
set_target_properties( ecrad_${PREC}_flux_helper_ranges PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module" )

if( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
    target_compile_options( ecrad_${PREC}_flux_helper_ranges PRIVATE "-mp=gpu" )
    target_link_options( ecrad_${PREC}_flux_helper_ranges PRIVATE "-mp=gpu" )
endif()
