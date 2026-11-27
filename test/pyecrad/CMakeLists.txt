# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Tests of the pyecrad Python package, using the libecrad4py library of
# each precision directly from the build tree (no installation):
#  - pytest tests of the Python interface (test_pyecrad.py)
#  - runs of the Python offline driver ("python -m pyecrad") compared
#    with runs of the Fortran offline driver on the same inputs, if
#    BITIDENTITY_TESTING is enabled (nccmp)

include( ecrad_python_macros )

# Optionally use the ARTIFACTS_DIR as wheelhouse, if provided
if( DEFINED ARTIFACTS_DIR )
    set( WHEELS_DIR_OPTION WHEELS_DIR "${ARTIFACTS_DIR}" )
else()
    set( WHEELS_DIR_OPTION "" )
endif()

# Python environment with the dependencies of pyecrad
ecrad_setup_python_venv( VENV_NAME pyecrad_venv )
ecrad_install_python_package(
    REQUIREMENT_SPEC -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt
    ${WHEELS_DIR_OPTION}
)
set( PYECRAD_PYTHON "${Python3_EXECUTABLE}" )

# Create symbolic link to the data used by the Fortran offline driver
file( CREATE_LINK "${PROJECT_SOURCE_DIR}/data" "${CMAKE_CURRENT_BINARY_DIR}/data" SYMBOLIC )

set( TEST_PYECRAD_IFS_INPUT "${PROJECT_SOURCE_DIR}/test/ifs/ecrad_meridian.nc" )
set( TEST_PYECRAD_IFS_CONFIG "${PROJECT_SOURCE_DIR}/test/ifs/configCY50R1.nam" )
set( TEST_PYECRAD_ECCKD_CONFIG "${PROJECT_SOURCE_DIR}/test/ifs/configCY50R1_ecckd.nam" )
set( TEST_PYECRAD_PRACTICAL_INPUT "${PROJECT_SOURCE_DIR}/practical/era5slice.nc" )
set( TEST_PYECRAD_PRACTICAL_CONFIG "${PROJECT_SOURCE_DIR}/practical/config.nam" )

function( add_pyecrad_test )
    # Run the Fortran (ecrad_${PREC}) and Python (python -m pyecrad)
    # offline drivers on INPUT_FILE with CONFIG_FILE modified by
    # NAMELIST_OPTIONS, and compare their outputs

    set( options "" )
    set( oneValueArgs NAME PREC CONFIG_FILE INPUT_FILE LW_THRESHOLD SW_THRESHOLD )
    set( multiValueArgs NAMELIST_OPTIONS )
    cmake_parse_arguments( _PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    set( _NAME pyecrad_${_PAR_PREC}_${_PAR_NAME} )

    # Outputs are written in double precision so that they can be
    # compared without loss of precision
    set( CONFIG_NAM "${CMAKE_CURRENT_BINARY_DIR}/config_${_NAME}.nam" )
    add_custom_command(
        OUTPUT ${CONFIG_NAM}
        COMMAND
            change_namelist.sh ${_PAR_CONFIG_FILE} ${CONFIG_NAM}
                "directory_name=\\\"data\\\"" do_write_double_precision=true
                ${_PAR_NAMELIST_OPTIONS}
        DEPENDS ${_PAR_CONFIG_FILE}
    )
    add_custom_target( config_${_NAME} ALL DEPENDS ${CONFIG_NAM} )

    ecbuild_add_test(
        TARGET ecrad_${_NAME}
        COMMAND ecrad_${_PAR_PREC}
        ARGS ${CONFIG_NAM} ${_PAR_INPUT_FILE} output_ecrad_${_NAME}.nc
        DEPENDS config_${_NAME}
    )

    ecbuild_add_test(
        TARGET ${_NAME}
        COMMAND ${PYECRAD_PYTHON}
        ARGS -m pyecrad ${CONFIG_NAM} ${_PAR_INPUT_FILE} output_${_NAME}.nc
        ENVIRONMENT ${PYECRAD_TEST_ENVIRONMENT}
        DEPENDS config_${_NAME} ecrad4py.${_PAR_PREC}
    )

    if( TARGET nccmp )
        ecbuild_add_test(
            TARGET nccmp_${_NAME}
            COMMAND nccmp
            ARGS
                --longwave-threshold=${_PAR_LW_THRESHOLD}
                --shortwave-threshold=${_PAR_SW_THRESHOLD}
                output_ecrad_${_NAME}.nc
                output_${_NAME}.nc
            TEST_DEPENDS ecrad_${_NAME} ${_NAME}
        )
    endif()

endfunction( add_pyecrad_test )

# Intel compilers call scalar math functions by their standard names
# (exp, log, pow...) and link the Intel math library (libimf) before the
# system libm: the Fortran offline driver uses the libimf versions. In the
# Python process, the system libm is loaded first by the interpreter and
# its versions are used instead, which changes the last bit of some
# results. Preloading libimf reproduces the symbol resolution of the
# Fortran driver so that both drivers can be compared bit for bit.
set( PYECRAD_PRELOAD_ENVIRONMENT "" )
if( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
    find_library( PYECRAD_INTEL_IMF_LIBRARY
        NAMES imf
        HINTS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}
    )
    if( PYECRAD_INTEL_IMF_LIBRARY )
        set( PYECRAD_PRELOAD_ENVIRONMENT "LD_PRELOAD=${PYECRAD_INTEL_IMF_LIBRARY}" )
    else()
        ecbuild_warn( "Intel math library (libimf) not found: the pyecrad results "
                      "may differ slightly from those of the Fortran offline driver" )
    endif()
endif()

foreach( PREC IN LISTS _PRECISIONS )

    # The Python package is taken from the source tree, and the library
    # and data files from the build tree
    set( PYECRAD_TEST_ENVIRONMENT
        "PYTHONPATH=${PROJECT_SOURCE_DIR}"
        "PYECRAD_LIBRARY=$<TARGET_FILE:ecrad4py.${PREC}>"
        "PYECRAD_DATA=${CMAKE_CURRENT_BINARY_DIR}/data"
        ${PYECRAD_PRELOAD_ENVIRONMENT}
    )

    # Tests of the Python interface
    ecbuild_add_test(
        TARGET pyecrad_${PREC}_pytest
        COMMAND ${PYECRAD_PYTHON}
        ARGS -m pytest -v -p no:cacheprovider ${CMAKE_CURRENT_SOURCE_DIR}/test_pyecrad.py
        ENVIRONMENT ${PYECRAD_TEST_ENVIRONMENT}
        DEPENDS ecrad4py.${PREC}
    )

    # Both drivers pass the same inputs to the radiation scheme, so the
    # results are expected to be identical
    set( LW_THRESHOLD 0.0 )
    set( SW_THRESHOLD 0.0 )

    # IFS CY50R1 configuration: McICA solver, RRTMG gas optics
    add_pyecrad_test(
        NAME default PREC ${PREC}
        CONFIG_FILE ${TEST_PYECRAD_IFS_CONFIG}
        INPUT_FILE ${TEST_PYECRAD_IFS_INPUT}
        LW_THRESHOLD ${LW_THRESHOLD} SW_THRESHOLD ${SW_THRESHOLD}
    )

    # Tripleclouds solver
    add_pyecrad_test(
        NAME tripleclouds PREC ${PREC}
        CONFIG_FILE ${TEST_PYECRAD_IFS_CONFIG}
        INPUT_FILE ${TEST_PYECRAD_IFS_INPUT}
        NAMELIST_OPTIONS
            "sw_solver_name=\\\"Tripleclouds\\\""
            "lw_solver_name=\\\"Tripleclouds\\\""
        LW_THRESHOLD ${LW_THRESHOLD} SW_THRESHOLD ${SW_THRESHOLD}
    )

    # SPARTACUS solver: the inverse cloud effective size is read from the
    # input file by both drivers (the cloud separation scale options of
    # the Fortran driver are not supported by pyecrad)
    add_pyecrad_test(
        NAME spartacus PREC ${PREC}
        CONFIG_FILE ${TEST_PYECRAD_IFS_CONFIG}
        INPUT_FILE ${TEST_PYECRAD_IFS_INPUT}
        NAMELIST_OPTIONS
            "sw_solver_name=\\\"SPARTACUS\\\""
            "lw_solver_name=\\\"SPARTACUS\\\""
            "do_3d_effects=true"
            "do_sw_delta_scaling_with_gases=false"
            "cloud_separation_scale_toa=-1.0"
            "cloud_separation_scale_surface=-1.0"
        LW_THRESHOLD ${LW_THRESHOLD} SW_THRESHOLD ${SW_THRESHOLD}
    )

    # No clouds and no aerosols
    add_pyecrad_test(
        NAME cloudless PREC ${PREC}
        CONFIG_FILE ${TEST_PYECRAD_IFS_CONFIG}
        INPUT_FILE ${TEST_PYECRAD_IFS_INPUT}
        NAMELIST_OPTIONS
            "use_aerosols=false"
            "sw_solver_name=\\\"Cloudless\\\""
            "lw_solver_name=\\\"Cloudless\\\""
        LW_THRESHOLD ${LW_THRESHOLD} SW_THRESHOLD ${SW_THRESHOLD}
    )

    # ecCKD gas optics
    add_pyecrad_test(
        NAME ecckd PREC ${PREC}
        CONFIG_FILE ${TEST_PYECRAD_ECCKD_CONFIG}
        INPUT_FILE ${TEST_PYECRAD_IFS_INPUT}
        LW_THRESHOLD ${LW_THRESHOLD} SW_THRESHOLD ${SW_THRESHOLD}
    )

    # Practical exercise: 1D surface properties, well-mixed oxygen and
    # overrides from the radiation_driver namelist. The overlap
    # decorrelation length scaling of 1 is disabled: the Fortran driver
    # computes overlap_param**1, which is not exact with approximate
    # math functions (e.g. Intel compilers with the fast floating-point
    # model), whereas numpy returns overlap_param.
    add_pyecrad_test(
        NAME practical PREC ${PREC}
        CONFIG_FILE ${TEST_PYECRAD_PRACTICAL_CONFIG}
        INPUT_FILE ${TEST_PYECRAD_PRACTICAL_INPUT}
        NAMELIST_OPTIONS
            "overlap_decorr_length_scaling=-1.0"
        LW_THRESHOLD ${LW_THRESHOLD} SW_THRESHOLD ${SW_THRESHOLD}
    )

endforeach()
