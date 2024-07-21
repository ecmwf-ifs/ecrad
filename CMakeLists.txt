# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

cmake_minimum_required( VERSION 3.20 FATAL_ERROR )

# ecrad uses ecbuild (https://github.com/ecmwf/ecbuild) to enhance the CMake build system.
# In an ifs-bundle build, ecbuild is already available as a shared dependency for multiple
# packages. Only if it cannot be found, ecbuild is downloaded automatically from the Github
# repository during the configure stage.

# For that, we first try to find ecbuild (QUIETly) ...
set( MINIMUM_ECBUILD_VERSION 3.8.5 )
find_package( ecbuild ${MINIMUM_ECBUILD_VERSION} HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild QUIET )

# ... or, if it cannot be found, download it via FetchContent ...
if( NOT ecbuild_FOUND )
    include(FetchContent)
    FetchContent_Declare(
        ecbuild
        GIT_REPOSITORY https://github.com/ecmwf/ecbuild.git
        GIT_TAG        ${MINIMUM_ECBUILD_VERSION}
    )
    FetchContent_MakeAvailable(ecbuild)
endif()

# ... before doing a final find_package to make it available.
find_package( ecbuild ${MINIMUM_ECBUILD_VERSION} REQUIRED )

# We set the <PNAME>_ECBUILD_COMPILE_FLAGS variable to an empty file. This ensures that we can
# specify a default set of compilation flags as well as per-file overwrites as part of the CMake config.
set( ECRAD_ECBUILD_COMPILE_FLAGS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/ecrad_init_flags.cmake CACHE FILEPATH "")

project( ecrad LANGUAGES Fortran VERSION 1.6.1 )
ecbuild_enable_fortran( REQUIRED )
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

set( ECRAD_COMPILE_DEFINITIONS "" )

# Dependencies
find_package( NetCDF COMPONENTS Fortran REQUIRED )

# Build options
ecbuild_add_option(
    FEATURE FIAT
    DEFAULT ON
    DESCRIPTION "Use PARKIND1, MPL, ABOR1, DR_HOOK from fiat"
    REQUIRED_PACKAGES "fiat"
)

ecbuild_add_option(
    FEATURE OMP
    DEFAULT ON
    DESCRIPTION "Use OpenMP for shared memory parallelisation in the driver programs"
    REQUIRED_PACKAGES "OpenMP COMPONENTS Fortran"
)

set( _PRECISIONS "" )
ecbuild_add_option(
    FEATURE DOUBLE_PRECISION
    DEFAULT ON
    DESCRIPTION "Build ecrad in double precision configuration"
)
if( HAVE_DOUBLE_PRECISION )
    list( APPEND _PRECISIONS "dp" )
endif()

ecbuild_add_option(
    FEATURE SINGLE_PRECISION
    DEFAULT OFF
    DESCRIPTION "Build ecrad in single precision configuration"
)
if( HAVE_SINGLE_PRECISION )
    list( APPEND _PRECISIONS "sp" )
endif()

ecbuild_add_option(
    FEATURE SHARED_LIBS
    DEFAULT ON
    DESCRIPTION "Build using shared libraries"
)
set( BUILD_SHARED_LIBS ${HAVE_SHARED_LIBS} )
if( HAVE_SHARED_LIBS )
    set( CMAKE_POSITION_INDEPENDENT_CODE TRUE )
endif()

ecbuild_add_option(
    FEATURE PRINT_ENTRAPMENT_DATA
    DEFAULT OFF
    DESCRIPTION "Write data to fort.101 and fort.102 in the SPARTACUS shortwave solver"
)
if( HAVE_PRINT_ENTRAPMENT_DATA )
    list( APPEND ECRAD_COMPILE_DEFINITIONS PRINT_ENTRAPMENT_DATA )
endif()

include(ecrad_compile_flags)

# Dummy drhook module
if( NOT HAVE_FIAT )
    add_subdirectory(drhook)
endif()

# Include subdirectories, once per precision, and using precision-specific
# binary directories
foreach( PREC IN LISTS _PRECISIONS )
    # Specify precision-specific Fortran module directory
    set( CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/radiation.${PREC}/module" )
    file( MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} )

    # Create an interface library to propagate external dependencies (OpenMP, FIAT),
    # or precision-specific compilation options to all targets
    ecbuild_add_library(
        TARGET ecrad_base.${PREC}
        TYPE INTERFACE
        PUBLIC_DEFINITIONS
            NC_NETCDF4
            ${ECRAD_COMPILE_DEFINITIONS}
            $<$<STREQUAL:"${PREC}","sp">:PARKIND1_SINGLE>
            $<${HAVE_FIAT}:HAVE_FIAT>
        PUBLIC_LIBS
            $<${HAVE_FIAT}:fiat parkind_${PREC}>
            $<$<NOT:${HAVE_FIAT}>:drhook_dummy>
    )

    add_subdirectory(ifsaux ifsaux.${PREC})
    add_subdirectory(utilities utilities.${PREC})
    add_subdirectory(ifsrrtm ifsrrtm.${PREC})
    add_subdirectory(radiation radiation.${PREC})
    add_subdirectory(ifs ifs.${PREC})
    add_subdirectory(driver driver.${PREC})
endforeach()

if( HAVE_TESTS )
    add_subdirectory(test)
endif()

ecbuild_install_project( NAME ${PROJECT_NAME} )
ecbuild_print_summary()
