# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Add change_namelist.sh as imported executable for use
# in setting up tests
add_executable( change_namelist.sh IMPORTED GLOBAL )
set_property(
    TARGET change_namelist.sh
    PROPERTY IMPORTED_LOCATION
        "${CMAKE_CURRENT_SOURCE_DIR}/change_namelist.sh"
)

# Set-up Python environment to enable the nccmp.py script
if( HAVE_BITIDENTITY_TESTING )
    include( ecrad_python_macros )

    # Optionally use the ARTIFACTS_DIR as wheelhouse, if provided
    if( DEFINED ARTIFACTS_DIR )
        set( WHEELS_DIR_OPTION WHEELS_DIR "${ARTIFACTS_DIR}" )
    else()
        set( WHEELS_DIR_OPTION "" )
    endif()

    # Install the requirements for nccmp.py
    ecrad_setup_python_venv( VENV_NAME nccmp_venv )
    ecrad_install_python_package(
        REQUIREMENT_SPEC -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt
        ${WHEELS_DIR_OPTION}
    )

    # Install a wrapper script to the build directory and import it as an executable
    configure_file(
        nccmp.in nccmp
        FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
        @ONLY
    )
    add_executable( nccmp IMPORTED GLOBAL )
    set_property( TARGET nccmp PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/nccmp )
endif()
