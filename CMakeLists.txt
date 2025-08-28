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
set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
include( ecrad_macros )


# ... or, if it cannot be found, download it via FetchContent ...
if( NOT ecbuild_FOUND )
    include(FetchContent)
    set( ECBUILD_ENABLE_TESTS OFF CACHE BOOL "Disable tests for ecbuild" )
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

project( ecrad LANGUAGES Fortran )
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

if( HAVE_FIAT AND fiat_HAVE_MPI )
    list( APPEND ECRAD_COMPILE_DEFINITIONS EASY_NETCDF_READ_MPI )
endif()

ecbuild_add_option(
    FEATURE OMP
    DEFAULT OFF
    DESCRIPTION "Use OpenMP for device offload or for shared memory parallelisation"
    REQUIRED_PACKAGES "OpenMP COMPONENTS Fortran"
)

if( ${CMAKE_VERSION} VERSION_LESS "3.25" )
  if ( ECRAD_ENABLE_ACC OR (NOT DEFINED ECRAD_ENABLE_ACC AND ENABLE_ACC) )
    # See https://gitlab.kitware.com/cmake/cmake/-/issues/23691, fixed in CMake 3.25
    # (TL;DR: FindOpenACC sets OpenACC_<LANG>_FOUND correctly but does not set
    #  OpenACC_FOUND unless all three C, CXX, and Fortran have been found - even if
    #  only one language has been requested via COMPONENTS)
    find_package( OpenACC COMPONENTS Fortran )
    if( OpenACC_Fortran_FOUND )
      set( OpenACC_FOUND ON )
    endif()
  endif()
endif()
ecbuild_add_option(
    FEATURE ACC
    DEFAULT ON
    DESCRIPTION "Use OpenACC to offload computations to GPU"
    REQUIRED_PACKAGES "OpenACC COMPONENTS Fortran"
    CONDITION CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC"
)
if( HAVE_ACC )
    if( DEFINED BUILD_SHARED_LIBS AND BUILD_SHARED_LIBS )
        ecbuild_warn( "OpenACC enabled: Switching to static linking" )
    endif()
    set( BUILD_SHARED_LIBS OFF CACHE BOOL "Build using shared libraries" FORCE )
endif()

ecbuild_add_option( FEATURE GPU
                    DEFAULT OFF
                    DESCRIPTION "Compile GPU version of ECRad (Requires OpenACC or sufficient OpenMP offloading support)"
                    CONDITION (HAVE_ACC OR HAVE_OMP) )

if( HAVE_GPU )
  if( HAVE_ACC )
    set( GPU_OFFLOAD "ACC" )
        ecbuild_warn( "OpenACC enabled" )
  elseif( HAVE_OMP )
    set( GPU_OFFLOAD "OMP" )
        ecbuild_warn( "OpenMP offload enabled" )
  else()
    ecbuild_error("Could not enable GPU as OMP or ACC were not enabled")
  endif()
endif()


if( CMAKE_VERSION VERSION_LESS 3.25 )
    set( NVTX_TARGET CUDA::nvToolsExt )
else()
    set( NVTX_TARGET CUDA::nvtx3 )
endif()
ecbuild_add_option(
    FEATURE NVTX
    DEFAULT ON
    DESCRIPTION "Enable NVTX markers in code"
    REQUIRED_PACKAGES "CUDAToolkit COMPONENTS ${NVTX_TARGET}"
    CONDITION HAVE_ACC
)
if( HAVE_NVTX )
    list( APPEND NVTX_TARGET nvhpcwrapnvtx )
    list( APPEND ECRAD_COMPILE_DEFINITIONS HAVE_NVTX )
endif()

if (HAVE_GPU)
  set(HAVE_HIP 0)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24")
    ecrad_find_hip() # sets "HAVE_HIP"
  else()
    message("Cannot enable for HIP language, require CMake >= 3.24")
  endif()
endif()

ecbuild_add_option(
  FEATURE ROCTX
  DEFAULT OFF
  DESCRIPTION "Enable Roctx markers in code"
  CONDITION HAVE_HIP
)

if( HAVE_ROCTX )
  list( APPEND ECRAD_COMPILE_DEFINITIONS HAVE_ROCTX )
endif()

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
    FEATURE PRINT_ENTRAPMENT_DATA
    DEFAULT OFF
    DESCRIPTION "Write data to fort.101 and fort.102 in the SPARTACUS shortwave solver"
)
if( HAVE_PRINT_ENTRAPMENT_DATA )
    list( APPEND ECRAD_COMPILE_DEFINITIONS PRINT_ENTRAPMENT_DATA )
endif()

ecbuild_add_option(
    FEATURE BITIDENTITY_TESTING
    DEFAULT OFF
    DESCRIPTION "Prescribe eff. radii, seed and cloud_fraction in blocked drivers for bit-identical results and validate with a Python script"
)
if( HAVE_BITIDENTITY_TESTING )
    list( APPEND ECRAD_COMPILE_DEFINITIONS BITIDENTITY_TESTING )
endif()

find_program( CDO NAMES cdo )
if( "${CDO}" STREQUAL "CDO-NOTFOUND" )
    ecbuild_info( "Could NOT find cdo" )
else()
    ecbuild_info( "Found program cdo [${CDO}]" )
    set( HAVE_CDO ON CACHE BOOL "cdo binary" )
endif()

include(ecrad_compile_flags)

# Dummy drhook module (precision-independent)
if( NOT HAVE_FIAT )
    add_subdirectory(drhook)
endif()

# Include subdirectories, once per precision, and using precision-specific
# binary directories
foreach( PREC IN LISTS _PRECISIONS )
    # Specify precision-specific Fortran module directory
    set( CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/radiation.${PREC}/module" )
    file( MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} )

    # Define precision-specific include directory for installing files
    set( ecrad_${PREC}_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/ecrad.${PREC}" )

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

# Include test directories only if tests are enabled (which is the default
# and can be disabled by setting ENABLE_TESTS=OFF)
if( HAVE_TESTS )
    add_subdirectory(test)
endif()

ecbuild_install_project( NAME ${PROJECT_NAME} )
ecbuild_print_summary()
