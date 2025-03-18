# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

set( TEST_IFS_INPUT "${CMAKE_CURRENT_SOURCE_DIR}/ecrad_meridian.nc" )
set( TEST_IFS_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/configCY49R1.nam" )
set( TEST_ECCKD_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/configCY49R1_ecckd.nam" )

# Create symbolic link to test data
file( CREATE_LINK "${PROJECT_SOURCE_DIR}/data" "${CMAKE_CURRENT_BINARY_DIR}/data" SYMBOLIC )

function( add_ecrad_nccmp_test )
    # Utility function to validate output from a test against a reference file

    set( options ENABLED XFAIL_VALIDATION )
    set( oneValueArgs NAME REFERENCE_NAME REFERENCE_FILE TEST_NAME LW_THRESHOLD SW_THRESHOLD )
    set( multiValueArgs )
    cmake_parse_arguments( _P "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if( TARGET nccmp AND "${PREC}" STREQUAL "dp" )

        if( _P_REFERENCE_FILE )
            if( _P_REFERENCE_NAME )
                ecbuild_error( "Options REFERENCE_NAME and REFERENCE_FILE are mutually exclusive" )
            endif()
        else()
            set( _P_REFERENCE_FILE output_${_P_REFERENCE_NAME}.nc )
        endif()

        set( DEPS ${_P_TEST_NAME} )
        if( _P_REFERENCE_NAME )
            list( APPEND DEPS ${_P_REFERENCE_NAME} )
        endif()

        ecbuild_add_test(
            TARGET nccmp_${_P_TEST_NAME}
            COMMAND nccmp
            ARGS
                --longwave-threshold=${_P_LW_THRESHOLD}
                --shortwave-threshold=${_P_SW_THRESHOLD}
                ${_P_REFERENCE_FILE}
                output_${_P_TEST_NAME}.nc
            TEST_DEPENDS ${DEPS}
        )

        if( ${_P_XFAIL_VALIDATION} )
            set_property( TEST nccmp_${_P_TEST_NAME} PROPERTY WILL_FAIL true )
        endif()

    endif()

endfunction( add_ecrad_nccmp_test )

function( add_ecrad_ifs_test )
    # Utility function to add tests for a specified configuration,
    # taking CONFIG_FILE as a template and modifying entries
    # using the provided NAMELIST_OPTIONS.
    #
    # If CONFIG FILE is not provided, it defaults to ${TEST_IFS_CONFIG}.
    #
    # This creates four tests for each configuration:
    #  1. Using ecrad.${PREC}
    #  2. Using ecrad.${PREC} with net flux outputs
    #  3. Using ecrad_ifs.${PREC} with net flux outputs
    #  4. Using ecrad_ifs_blocked.${PREC} with net flux outputs

    set( options ENABLED XFAIL_VALIDATION )
    set( oneValueArgs NAME CONFIG_FILE LW_THRESHOLD SW_THRESHOLD )
    set( multiValueArgs NAMELIST_OPTIONS )
    cmake_parse_arguments( _PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if( ${_PAR_ENABLED} )
        # Setup namelist file for the test, either by overwriting options
        # or by copying the original file
        if( NOT _PAR_CONFIG_FILE )
            set( _PAR_CONFIG_FILE ${TEST_IFS_CONFIG} )
        endif()
        set( CONFIG_NAM "${CMAKE_CURRENT_BINARY_DIR}/config_${_PAR_NAME}.nam" )
        list( APPEND _PAR_NAMELIST_OPTIONS "directory_name=\\\"data\\\"" )
        add_custom_command(
            OUTPUT ${CONFIG_NAM}
            COMMAND
                change_namelist.sh ${_PAR_CONFIG_FILE} ${CONFIG_NAM} ${_PAR_NAMELIST_OPTIONS}
            DEPENDS ${_PAR_CONFIG_FILE}
        )
        add_custom_target( config_${_PAR_NAME}_${PREC} DEPENDS ${CONFIG_NAM} )

        # Create second namelist file with net fluxes enabled
        set( CONFIG_NET_NAM "${CMAKE_CURRENT_BINARY_DIR}/config_${_PAR_NAME}_net.nam" )
        add_custom_command(
            OUTPUT ${CONFIG_NET_NAM}
            COMMAND
                change_namelist.sh ${CONFIG_NAM} ${CONFIG_NET_NAM}
                    do_save_net_fluxes=true do_write_double_precision=true
            DEPENDS config_${_PAR_NAME}_${PREC}
        )
        add_custom_target( config_${_PAR_NAME}_${PREC}_net DEPENDS ${CONFIG_NET_NAM} )

        ecbuild_add_test(
            TARGET ecrad_${PREC}_${_PAR_NAME}
            COMMAND ecrad_${PREC}
            ARGS
                ${CONFIG_NAM}
                ${TEST_IFS_INPUT}
                output_ecrad_${PREC}_${_PAR_NAME}.nc
            DEPENDS config_${_PAR_NAME}_${PREC}
        )

        foreach( binary ecrad_${PREC} ecrad_ifs_${PREC} ecrad_ifs_blocked_${PREC} )

            ecbuild_add_test(
                TARGET ${binary}_${_PAR_NAME}_net
                COMMAND ${binary}
                ARGS
                    ${CONFIG_NET_NAM}
                    ${TEST_IFS_INPUT}
                    output_${binary}_${_PAR_NAME}_net.nc
                DEPENDS config_${_PAR_NAME}_${PREC}_net
            )

        endforeach()

        # When bit-identity testing is enabled, we can check the output against the
        # default standalone driver
        if( HAVE_BITIDENTITY_TESTING )

            if( NOT _PAR_LW_THRESHOLD )
                set( _PAR_LW_THRESHOLD 0.0 )
            endif()
            if( NOT _PAR_SW_THRESHOLD )
                set( _PAR_SW_THRESHOLD 0.0 )
            endif()

            if( _PAR_XFAIL_VALIDATION )
                set( NCCMP_OPTIONS XFAIL_VALIDATION )
            else()
                set( NCCMP_OPTIONS )
            endif()

            foreach( binary ecrad_ifs_${PREC} ecrad_ifs_blocked_${PREC} )

                if( TARGET ${binary} )
                    add_ecrad_nccmp_test(
                        REFERENCE_NAME ecrad_${PREC}_${_PAR_NAME}_net
                        TEST_NAME ${binary}_${_PAR_NAME}_net
                        LW_THRESHOLD ${_PAR_LW_THRESHOLD}
                        SW_THRESHOLD ${_PAR_SW_THRESHOLD}
                        ${NCCMP_OPTIONS}
                    )
                endif()

            endforeach()

        endif()

    endif()
endfunction( add_ecrad_ifs_test )

##############################################################################
# Tests of CY49R1 configuration with RRTMG gas optics model and
# various solvers, using ${TEST_IFS_CONFIG} config file

# IFS default: McICA solver with exponential-exponential overlap
add_ecrad_ifs_test(
    NAME default
    ENABLED ON
)

# When the nccmp utility is available, we can check the output against the reference file
# in the repository
add_ecrad_nccmp_test(
    REFERENCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/ecrad_meridian_default_out_REFERENCE.nc
    TEST_NAME ecrad_${PREC}_default
    LW_THRESHOLD 0.001
    SW_THRESHOLD 0.1
)

# Turn off aerosols
add_ecrad_ifs_test(
    NAME noaer
    NAMELIST_OPTIONS "use_aerosols=false"
    ENABLED ON
)

# Older exponential-exponential overlap
add_ecrad_ifs_test(
    NAME expexp
    NAMELIST_OPTIONS "overlap_scheme_name=\\\"Exp-Exp\\\""
    ENABLED ON
)

# Tripleclouds solver with exponential-random overlap
add_ecrad_ifs_test(
    NAME tripleclouds
	NAMELIST_OPTIONS
        "sw_solver_name=\\\"Tripleclouds\\\""
        "lw_solver_name=\\\"Tripleclouds\\\""
    ENABLED ON
)

# Longwave scattering; since 46R1 this is the default
add_ecrad_ifs_test(
    NAME lwscat
    NAMELIST_OPTIONS "do_lw_cloud_scattering=true"
    ENABLED ON
)

#  3D radiative transfer
add_ecrad_ifs_test(
    NAME spartacus
    NAMELIST_OPTIONS
		"sw_solver_name=\\\"SPARTACUS\\\""
        "lw_solver_name=\\\"SPARTACUS\\\""
		"do_3d_effects=true"
		"do_sw_delta_scaling_with_gases=false"
    ENABLED ON
    XFAIL_VALIDATION
)

# 3D radiative transfer using the older "maximum entrapment"
add_ecrad_ifs_test(
    NAME spartacus_maxentr
    NAMELIST_OPTIONS
		"sw_solver_name=\\\"SPARTACUS\\\""
        "lw_solver_name=\\\"SPARTACUS\\\""
		"do_3d_effects=true"
		"sw_entrapment_name=\\\"Maximum\\\""
		"do_sw_delta_scaling_with_gases=false"
    ENABLED ON
    XFAIL_VALIDATION
)

# "Cloudless" solver
add_ecrad_ifs_test(
    NAME cloudless
    NAMELIST_OPTIONS
		"use_aerosols=false"
		"sw_solver_name=\\\"Cloudless\\\""
        "lw_solver_name=\\\"Cloudless\\\""
    ENABLED ON
)

# Exponential-random overlap with "vectorizable" cloud generator
add_ecrad_ifs_test(
    NAME vec
    NAMELIST_OPTIONS "use_vectorizable_generator=true"
    ENABLED ON
)

##############################################################################
# The following targets use the ${TEST_ECCKD_CONFIG} configuration file

# ecCKD gas optics scheme (note that default solver is Tripleclouds)
add_ecrad_ifs_test(
    NAME ecckd_mcica
    CONFIG_FILE "${TEST_ECCKD_CONFIG}"
    NAMELIST_OPTIONS
        "sw_solver_name=\\\"McICA\\\""
        "lw_solver_name=\\\"McICA\\\""
    SW_THRESHOLD 1e-10
    LW_THRESHOLD 1e-10
    ENABLED ON
)

# ecCKD with Tripleclouds solver (default)
add_ecrad_ifs_test(
    NAME ecckd_tc
    CONFIG_FILE "${TEST_ECCKD_CONFIG}"
    SW_THRESHOLD 1e-10
    LW_THRESHOLD 1e-10
    ENABLED ON
)

# ecCKD with no aerosols
add_ecrad_ifs_test(
    NAME ecckd_noaer
    CONFIG_FILE "${TEST_ECCKD_CONFIG}"
    NAMELIST_OPTIONS "use_aerosols=false"
    SW_THRESHOLD 1e-10
    LW_THRESHOLD 1e-10
    ENABLED ON
)

# ecCKD gas optics with SPARTACUS solver (not currently correct)
add_ecrad_ifs_test(
    NAME ecckd_spartacus
    CONFIG_FILE "${TEST_ECCKD_CONFIG}"
    NAMELIST_OPTIONS
        "sw_solver_name=\\\"SPARTACUS\\\""
        "lw_solver_name=\\\"SPARTACUS\\\""
        "do_3d_effects=true"
    ENABLED ON
    XFAIL_VALIDATION
)
