# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

set( radiation_SOURCES
        radiation_adding_ica_lw.F90
        radiation_adding_ica_sw.F90
        radiation_aerosol_optics_data.F90
        radiation_aerosol_optics_description.F90
        radiation_aerosol_optics.F90
        radiation_aerosol.F90
        radiation_check.F90
        radiation_cloud_cover.F90
        radiation_cloud_generator.F90
        radiation_cloud_generator_acc.F90
        radiation_cloud_optics_data.F90
        radiation_cloud_optics.F90
        radiation_cloud.F90
        radiation_cloudless_lw.F90
        radiation_cloudless_sw.F90
        radiation_config.F90
        radiation_constants.F90
        radiation_delta_eddington.h
        radiation_ecckd_gas.F90
        radiation_ecckd_interface.F90
        radiation_ecckd.F90
        radiation_flux.F90
        radiation_gas_constants.F90
        radiation_gas.F90
        radiation_general_cloud_optics_data.F90
        radiation_general_cloud_optics.F90
        radiation_homogeneous_lw.F90
        radiation_homogeneous_sw.F90
        radiation_ice_optics_baran.F90
        radiation_ice_optics_baran2016.F90
        radiation_ice_optics_baran2017.F90
        radiation_ice_optics_fu.F90
        radiation_ice_optics_yi.F90
        radiation_ifs_rrtm.F90
        radiation_interface.F90
        radiation_liquid_optics_slingo.F90
        radiation_liquid_optics_socrates.F90
        radiation_lw_derivatives.F90
        radiation_matrix.F90
        radiation_mcica_acc_lw.F90
        radiation_mcica_acc_sw.F90
        radiation_mcica_lw.F90
        radiation_mcica_sw.F90
        radiation_monochromatic.F90
        radiation_optical_depth_scaling.h
        radiation_overlap.F90
        radiation_pdf_sampler.F90
        radiation_random_numbers.F90
        radiation_regions.F90
        radiation_save.F90
        radiation_single_level.F90
        radiation_spartacus_lw.F90
        radiation_spartacus_sw.F90
        radiation_spectral_definition.F90
        radiation_thermodynamics.F90
        radiation_tripleclouds_lw.F90
        radiation_tripleclouds_sw.F90
        radiation_two_stream.F90
)

ecbuild_add_library(
    TARGET
        ecrad.${PREC}
    SOURCES
        ${radiation_SOURCES}
    PUBLIC_INCLUDES
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
        "$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>"
        "$<INSTALL_INTERFACE:${ecrad_${PREC}_INSTALL_INCLUDEDIR}>"
    PRIVATE_LIBS
        ecrad_base.${PREC}
        ecrad_ifsaux.${PREC}
        ecrad_utilities.${PREC}
        ifsrrtm.${PREC}
        $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
    PUBLIC_LIBS
        $<${HAVE_ACC}:OpenACC::OpenACC_Fortran>
)

add_library( ecrad::ecrad.${PREC} ALIAS ecrad.${PREC} )
install(
    FILES
        "${CMAKE_Fortran_MODULE_DIRECTORY}/radiation_aerosol.mod"
        "${CMAKE_Fortran_MODULE_DIRECTORY}/radiation_aerosol_optics.mod"
        "${CMAKE_Fortran_MODULE_DIRECTORY}/radiation_aerosol_optics_description.mod"
        "${CMAKE_Fortran_MODULE_DIRECTORY}/radiation_aerosol_optics_data.mod"
        "${CMAKE_Fortran_MODULE_DIRECTORY}/radiation_cloud.mod"
        "${CMAKE_Fortran_MODULE_DIRECTORY}/radiation_cloud_cover.mod"
        "${CMAKE_Fortran_MODULE_DIRECTORY}/radiation_config.mod"
        "${CMAKE_Fortran_MODULE_DIRECTORY}/radiation_flux.mod"
        "${CMAKE_Fortran_MODULE_DIRECTORY}/radiation_gas.mod"
        "${CMAKE_Fortran_MODULE_DIRECTORY}/radiation_gas_constants.mod"
        "${CMAKE_Fortran_MODULE_DIRECTORY}/radiation_interface.mod"
        "${CMAKE_Fortran_MODULE_DIRECTORY}/radiation_save.mod"
        "${CMAKE_Fortran_MODULE_DIRECTORY}/radiation_single_level.mod"
        "${CMAKE_Fortran_MODULE_DIRECTORY}/radiation_thermodynamics.mod"
    DESTINATION ${ecrad_${PREC}_INSTALL_INCLUDEDIR}
)

# Compile flag overwrites
if( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
    set_source_files_properties( radiation_mcica_sw.F90 PROPERTIES COMPILE_OPTIONS "-vecabi=cmdtarget")
    set_source_files_properties( radiation_cloud_generator.F90 PROPERTIES COMPILE_OPTIONS "-vecabi=cmdtarget")
endif()
