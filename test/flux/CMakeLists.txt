# (C) Copyright 2026- ECMWF.
# Licensed under the Apache Licence Version 2.0; see COPYING.

ecbuild_add_test(
    TARGET ecrad_${PREC}_flux_helper_ranges
    SOURCES flux_helper_ranges.F90
    LIBS ecrad.${PREC} ecrad_base.${PREC} OpenMP::OpenMP_Fortran
)

if( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
    target_compile_options( ecrad_${PREC}_flux_helper_ranges PRIVATE "-mp=gpu" )
    target_link_options( ecrad_${PREC}_flux_helper_ranges PRIVATE "-mp=gpu" )
endif()
