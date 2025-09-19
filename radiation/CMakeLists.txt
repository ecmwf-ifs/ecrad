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
        radiation_liquid_optics_jahangir.F90
        radiation_liquid_optics_nielsen.F90
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

if( BUILD_SHARED_LIBS )
    if( HAVE_GPU )
        ecbuild_warn( "GPU support enabled: Forcing static linking for ecrad.${PREC}" )
        set( LIBRARY_TYPE STATIC )
    else()
        set( LIBRARY_TYPE SHARED )
    endif()
else()
    set( LIBRARY_TYPE STATIC )
endif()

ecbuild_add_library(
    TARGET
        ecrad.${PREC}
    TYPE
        ${LIBRARY_TYPE}
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
        $<${HAVE_NVTX}:${NVTX_TARGET}>
        $<${HAVE_ROCTX}:${ROCTX_LIBRARIES}>
    PRIVATE_DEFINITIONS ${GPU_OFFLOAD}GPU
)

# To ensure device code is generated and all link time dependencies are
# available with OpenMP target offloading, force the `-mp=gpu` flag for NVHPC
if( HAVE_GPU AND HAVE_OMP AND GPU_OFFLOAD STREQUAL "OMP" )
    if( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
        target_compile_options( ecrad.${PREC} PRIVATE "-mp=gpu" )
        target_link_options( ecrad.${PREC} PRIVATE "-mp=gpu" )
    endif()
endif()

# Improve performance via function inlining for NVHPC
# This requires building an inline lib for Fortran, see
# https://docs.nvidia.com/hpc-sdk/compilers/hpc-compilers-user-guide/index.html#fn-inline-use
set( inlib_NAME ecrad.${PREC}.il )
if( HAVE_GPU AND CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC" )
    set( HAVE_INLIB ON )
else()
    set( HAVE_INLIB OFF )
endif()

if( HAVE_INLIB )

    # List of source files from which functions should be inlined
    set( inlib_SOURCES
            radiation_random_numbers.F90
            radiation_lw_derivatives.F90
            radiation_adding_ica_lw.F90
            radiation_adding_ica_sw.F90
            radiation_cloud_generator_acc.F90
            radiation_two_stream.F90
    )
    set( inlib_FUNCTIONS
            initialize_acc uniform_distribution_acc
            calc_lw_derivatives_ica2 modify_lw_derivatives_ica
            adding_ica_lw fast_adding_ica_lw calc_fluxes_no_scattering_lw
            adding_ica_sw
            cloud_generator_acc
            calc_ref_trans_sw calc_ref_trans_lw calc_reflectance_transmittance_sw
                calc_no_scattering_transmittance_lw calc_two_stream_gammas_sw
    )

    # Additional sources that are required to compile inlib
    list( APPEND inlib_SOURCES
            radiation_matrix.F90
    )

    # Build the list of consumers
    set( inlib_CONSUMERS ${radiation_SOURCES} )
    list( REMOVE_ITEM inlib_CONSUMERS ${inlib_SOURCES} )

    # Declare the inlib library as an object library, although the extraction
    # pass does not actually generate object files. This makes sure we end up
    # using the same compilation options as for the actual ecrad library.
    # The ecrad library does not ingest the inlib library as a LIBRARY dependency,
    # but instead we add an explicit inter-target dependency with add_dependencies
    # to ensure the inline library is being built before the ecrad library.
    ecbuild_add_library(
        TARGET ${inlib_NAME}
        TYPE OBJECT
        SOURCES
            ${inlib_SOURCES}
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

    # To ensure device code is generated and all link time dependencies are
    # available with OpenMP target offloading, force the `-mp=gpu` flag for NVHPC
    if( HAVE_GPU AND HAVE_OMP AND GPU_OFFLOAD STREQUAL "OMP" )
        if( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
            target_compile_options( ${inlib_NAME} PUBLIC "-mp=gpu" )
            target_link_options( ${inlib_NAME} PUBLIC "-mp=gpu" )
        endif()
    endif()

    # Make sure the inline library is created first
    add_dependencies( ecrad.${PREC} ${inlib_NAME} )

    # The compile flags to add to the relevant source files
    list( JOIN inlib_FUNCTIONS "," inlib_FUNCTIONS_CSV )
    set( inlib_extract_FCFLAGS "-Mextract=lib:${inlib_NAME},reshape,${inlib_FUNCTIONS_CSV}" )
    set( inlib_inline_FCFLAGS "-Minline=lib:${inlib_NAME},reshape,${inlib_FUNCTIONS_CSV};-Minfo=inline" )

    # Add compile options to extract functions to inline
    target_compile_options( ${inlib_NAME} PRIVATE ${inlib_extract_FCFLAGS} )

    # Add compile options to ingest inlined functions to relevant source files
    set_source_files_properties( ${inlib_CONSUMERS} PROPERTIES COMPILE_OPTIONS "${inlib_inline_FCFLAGS}" )

endif()

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
elseif( CMAKE_FORTRAN_COMPILER_ID MATCHES "PGI|NVHPC" )
    if( HAVE_GPU )
        set_source_files_properties(
            radiation_mcica_acc_lw.F90
            radiation_mcica_acc_sw.F90
            radiation_two_stream.F90
            radiation_lw_derivatives.F90
            PROPERTIES COMPILE_OPTIONS
            "-gpu=maxregcount:96"
        )
    endif()
endif()
