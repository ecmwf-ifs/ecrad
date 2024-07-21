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

    set( options ENABLED )
    set( oneValueArgs NAME CONFIG_FILE )
    set( multiValueArgs NAMELIST_OPTIONS )
    cmake_parse_arguments( _PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if( ${_PAR_ENABLED} )
        # Setup namelist file for the test, either by overwriting options
        # or by copying the original file
        set( CONFIG_NAM "${CMAKE_CURRENT_BINARY_DIR}/config_${_PAR_NAME}.nam" )
        list( APPEND _PAR_NAMELIST_OPTIONS "directory_name=\"data\"" )
        add_custom_command(
            OUTPUT ${CONFIG_NAM}
            COMMAND
                change_namelist.sh ${TEST_IFS_CONFIG} ${CONFIG_NAM} ${_PAR_NAMELIST_OPTIONS}
            DEPENDS ${TEST_IFS_CONFIG}
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
            TARGET ecrad.${PREC}_${_PAR_NAME}
            COMMAND ecrad.${PREC}
            ARGS
                ${CONFIG_NAM}
                ${TEST_IFS_INPUT}
                output_${_PAR_NAME}.nc
            DEPENDS config_${_PAR_NAME}_${PREC}
        )

        foreach( binary ecrad.${PREC} ecrad_ifs.${PREC} ecrad_ifs_blocked.${PREC} )

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

# Turn off aerosols
add_ecrad_ifs_test(
    NAME noaer
    NAMELIST_OPTIONS "use_aerosols=false"
    ENABLED ON
)

# Older exponential-exponential overlap
add_ecrad_ifs_test(
    NAME expexp
    NAMELIST_OPTIONS "overlap_scheme_name=\"Exp-Exp\""
    ENABLED ON
)

# Tripleclouds solver with exponential-random overlap
add_ecrad_ifs_test(
    NAME tripleclouds
	NAMELIST_OPTIONS
        "sw_solver_name=\"Tripleclouds\""
        "lw_solver_name=\"Tripleclouds\""
    ENABLED ON
)

# Longwave scattering; since 46R1 this is the default
add_ecrad_ifs_test(
    NAME lwscat
    NAMELIST_OPTIONS "do_lw_cloud_scattering=\"true\""
    ENABLED ON
)

#  3D radiative transfer
add_ecrad_ifs_test(
    NAME spartacus
    NAMELIST_OPTIONS
		"sw_solver_name=\"SPARTACUS\""
        "lw_solver_name=\"SPARTACUS\""
		"do_3d_effects=\"true\""
		"do_sw_delta_scaling_with_gases=\"false\""
    ENABLED ON
)

# 3D radiative transfer using the older "maximum entrapment"
add_ecrad_ifs_test(
    NAME spartacus_maxentr
    NAMELIST_OPTIONS
		"sw_solver_name=\"SPARTACUS\""
        "lw_solver_name=\"SPARTACUS\""
		"do_3d_effects=\"true\""
		"sw_entrapment_name=\"Maximum\""
		"do_sw_delta_scaling_with_gases=\"false\""
    ENABLED ON
)

# "Cloudless" solver
add_ecrad_ifs_test(
    NAME cloudless
    NAMELIST_OPTIONS
		"use_aerosols=false"
		"sw_solver_name=\"Cloudless\""
        "lw_solver_name=\"Cloudless\""
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
        "sw_solver_name=\"McICA\""
        "lw_solver_name=\"McICA\""
    ENABLED ON
)

# ecCKD with Tripleclouds solver (default)
add_ecrad_ifs_test(
    NAME ecckd_tc
    CONFIG_FILE "${TEST_ECCKD_CONFIG}"
    ENABLED ON
)

# ecCKD with no aerosols
add_ecrad_ifs_test(
    NAME ecckd_noaer
    CONFIG_FILE "${TEST_ECCKD_CONFIG}"
    NAMELIST_OPTIONS "use_aerosols=false"
    ENABLED ON
)

# ecCKD gas optics with SPARTACUS solver (not currently correct)
add_ecrad_ifs_test(
    NAME ecckd_spartacus
    CONFIG_FILE "${TEST_ECCKD_CONFIG}"
    NAMELIST_OPTIONS
        "sw_solver_name=\"SPARTACUS\""
        "lw_solver_name=\"SPARTACUS\""
        "do_3d_effects=\"true\""
    ENABLED ON
)
