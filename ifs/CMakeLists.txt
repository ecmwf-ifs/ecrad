# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

set( ifs_SOURCES
        ice_effective_radius.F90
        liquid_effective_radius.F90
        radiation_scheme.F90
        radiation_setup.F90
        yoerdu.F90
        yomrip.F90
        yoephy.F90
        yoecld.F90
        yoe_spectral_planck.F90
        cloud_overlap_decorr_len.F90
        yoerad.F90
        ifs_blocking.F90
${CMAKE_CURRENT_BINARY_DIR}/radintg_zrgp_mod.F90
)

macro( ecrad_ifs_process_fypp fypp_file )
    # If FYPP has been enabled, regenerate the F90 source files otherwise simply link
    # or copy commited versions
    if( HAVE_FYPP )
        add_custom_command(
            OUTPUT
                ${CMAKE_CURRENT_BINARY_DIR}/${fypp_file}.F90
            COMMAND
                ${FYPP} -m os -m field_config -M ${CMAKE_CURRENT_SOURCE_DIR}
                ${CMAKE_CURRENT_SOURCE_DIR}/${fypp_file}.fypp
                ${CMAKE_CURRENT_BINARY_DIR}/${fypp_file}.F90
            DEPENDS
                ${CMAKE_CURRENT_SOURCE_DIR}/${fypp_file}.fypp
                ${CMAKE_CURRENT_SOURCE_DIR}/radiation_fields_config.yaml
        )
    else()
        file( CREATE_LINK
            ${CMAKE_CURRENT_SOURCE_DIR}/${fypp_file}.F90
            ${CMAKE_CURRENT_BINARY_DIR}/${fypp_file}.F90
            SYMBOLIC COPY_ON_ERROR
        )
    endif()
endmacro()

ecrad_ifs_process_fypp( radintg_zrgp_mod )

ecbuild_add_library(
    TARGET ifs.${PREC}
    TYPE OBJECT
    SOURCES ${ifs_SOURCES}
    PUBLIC_INCLUDES
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/module_ifs>"
    PRIVATE_LIBS
        ecrad.${PREC}
        ecrad_base.${PREC}
)
set_target_properties( ifs.${PREC}
    PROPERTIES
        Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module_ifs"
)
