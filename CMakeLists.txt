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

ecbuild_add_option(
    FEATURE PYECRAD
    DEFAULT OFF
    DESCRIPTION "Build and install the pyecrad Python wrapper in a virtual environment"
)

option( ENABLE_STATIC_LIBS
    "Build static libraries for ecrad and drhook_dummy" OFF
)

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

# ── pyecrad Python package ──────────────────────────────────────────────
#
# When enabled, sets up a Python virtual environment with pyecrad's
# dependencies and installs pyecrad as an editable package.  The ecrad4py
# shared library is staged into the source tree by a POST_BUILD command
# in driver/CMakeLists.txt and is picked up by pip install.
#
if( HAVE_PYECRAD )
    include( ecrad_python_macros )

    ecrad_setup_python_venv( VENV_NAME pyecrad_env )

    # Install pure-Python dependencies (the Fortran .so is already staged)
    # setuptools/wheel are needed for --no-build-isolation in the pip install step
    ecrad_install_python_package(
        REQUIREMENT_SPEC
            setuptools wheel
            "ctypesForFortran>=1.3.0, !=2.0.*, !=2.1.*"
            numpy f90nml netCDF4
    )

    # Collect the ecrad4py shared-library targets that were built
    set( _PYECRAD_DEPS "" )
    foreach( PREC IN LISTS _PRECISIONS )
        if( TARGET ecrad4py.${PREC} )
            list( APPEND _PYECRAD_DEPS ecrad4py.${PREC} )
        endif()
    endforeach()

    # Build-time: pip install pyecrad from the project root.
    # --no-build-isolation ensures pip uses the venv's build tools and
    # does not try to rebuild the .so (which CMake already handled).
    add_custom_target( pyecrad ALL
        COMMAND ${CMAKE_COMMAND} -E env
            "VIRTUAL_ENV=${CMAKE_CURRENT_BINARY_DIR}/pyecrad_env"
            "ECRAD_CMAKE_BUILD=1"
            "${Python3_VENV_BIN}/python" -m pip install
                --no-build-isolation -e "${CMAKE_SOURCE_DIR}"
        DEPENDS ${_PYECRAD_DEPS}
        COMMENT "Installing pyecrad in virtual environment"
    )

    # Export the Python executable for downstream targets (e.g. tests)
    set( PYECRAD_PYTHON "${Python3_EXECUTABLE}" CACHE INTERNAL "" )
endif()

# Include test directories only if tests are enabled (which is the default
# and can be disabled by setting ENABLE_TESTS=OFF)
if( HAVE_TESTS )
    add_subdirectory(test)
endif()

ecbuild_install_project( NAME ${PROJECT_NAME} )
ecbuild_print_summary()
