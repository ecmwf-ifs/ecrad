# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

ecbuild_add_library(
    TARGET driver_lib.${PREC}
    TYPE OBJECT
    SOURCES
        ecrad_driver_read_input.F90
        ecrad_driver_config.F90
        print_matrix_mod.F90
    PUBLIC_DEFINITIONS
        $<$<NOT:${HAVE_OMP}>:NO_OPENMP>
    PUBLIC_LIBS
        ecrad.${PREC}
        ecrad_base.${PREC}
        $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
)

ecbuild_add_executable(
    TARGET ecrad_${PREC}
    SOURCES
        ecrad_driver.F90
    LIBS
        driver_lib.${PREC}
)

ecbuild_add_executable(
    TARGET ecrad_ifs_${PREC}
    SOURCES
        ecrad_ifs_driver.F90
    LIBS
        ifs.${PREC}
        driver_lib.${PREC}
)

ecbuild_add_executable(
    TARGET ecrad_ifs_blocked_${PREC}
    SOURCES
        ecrad_ifs_driver_blocked.F90
        ifs_blocking.F90
    LIBS
        ifs.${PREC}
        driver_lib.${PREC}
)

if( ENABLE_STATIC_LIBS )
    set( CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/radiation.${PREC}/module_static" )
    add_library(ecrad4py.${PREC} SHARED
        ecrad4py.F90
        $<TARGET_OBJECTS:driver_lib.${PREC}>
        $<TARGET_OBJECTS:ifs.${PREC}>
    )
    target_compile_definitions(ecrad4py.${PREC} PRIVATE
        NC_NETCDF4
        ${ECRAD_COMPILE_DEFINITIONS}
        $<$<STREQUAL:"${PREC}","sp">:PARKIND1_SINGLE>
        $<${HAVE_FIAT}:HAVE_FIAT>
    )
    target_include_directories(ecrad4py.${PREC} PRIVATE
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
        "$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/radiation.${PREC}/module_static>"
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/driver>"
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/ifs>"
    )
    target_link_libraries(ecrad4py.${PREC} PRIVATE
        ecrad_static.${PREC}
        $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
    )
    if( NOT HAVE_FIAT )
        target_include_directories(ecrad4py.${PREC} PRIVATE
            "$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/drhook/static_mod>"
        )
        target_link_libraries(ecrad4py.${PREC} PRIVATE drhook_dummy_static)
    else()
        target_link_libraries(ecrad4py.${PREC} PRIVATE ecrad_base.${PREC})
    endif()
else()
    ecbuild_add_library(
        TARGET ecrad4py.${PREC}
        TYPE SHARED
        SOURCES ecrad4py.F90
        PRIVATE_LIBS
            driver_lib.${PREC}
            ifs.${PREC}
    )
endif()

# Set output name and stage the .so into the pyecrad source tree
if( "${PREC}" STREQUAL "dp" )
    set( _ECRAD4PY_SO "libecrad4py${CMAKE_SHARED_LIBRARY_SUFFIX}" )
    set_target_properties( ecrad4py.${PREC} PROPERTIES OUTPUT_NAME ecrad4py )
else()
    set( _ECRAD4PY_SO "libecrad4py_${PREC}${CMAKE_SHARED_LIBRARY_SUFFIX}" )
    set_target_properties( ecrad4py.${PREC} PROPERTIES OUTPUT_NAME ecrad4py_${PREC} )
endif()
add_custom_command( TARGET ecrad4py.${PREC} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:ecrad4py.${PREC}>
        ${PROJECT_SOURCE_DIR}/pyecrad/${_ECRAD4PY_SO}
    COMMENT "Staging ${_ECRAD4PY_SO} to pyecrad/"
)

# For ecrad_ifs_driver_blocked we have to disable bounds checking
# (which is enabled in Debug builds) because sequence association
# and values initialized to -9999 for unused arguments is used in
# the call to RADIATION_SCHEME.
if( CMAKE_BUILD_TYPE MATCHES "Debug" )
    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        set_source_files_properties(
            "ecrad_ifs_driver_blocked.F90"
            PROPERTIES COMPILE_FLAGS "-fcheck=no-bounds"
        )
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        set_source_files_properties(
            "ecrad_ifs_driver_blocked.F90"
            PROPERTIES COMPILE_FLAGS "-check nobounds"
        )
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC")
        set_source_files_properties(
            "ecrad_ifs_driver_blocked.F90"
            PROPERTIES COMPILE_FLAGS "-Mnobounds"
        )
    endif()
endif()

# Create a symlink for each driver without the precision suffix.
# The default is to link to the double precision version if it is built,
# and single precision otherwise.
if ( "${PREC}" STREQUAL "dp" OR NOT HAVE_DOUBLE_PRECISION )
    foreach( tgt ecrad ecrad_ifs ecrad_ifs_blocked )
        add_custom_command(
            TARGET ${tgt}_${PREC}
            POST_BUILD
            COMMAND
                "${CMAKE_COMMAND}" -E create_symlink ${tgt}_${PREC} ${tgt}
            BYPRODUCTS ${tgt}
            WORKING_DIRECTORY
                "${CMAKE_BINARY_DIR}/bin"
            COMMENT "Creating symbolic link from ${tgt} to ${tgt}_${PREC}"
        )
        install(
            FILES ${CMAKE_BINARY_DIR}/bin/${tgt}
            DESTINATION bin
        )
    endforeach()
endif()
