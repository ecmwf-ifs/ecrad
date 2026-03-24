# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Add here extra compile flags for specific files
# This file gets included in the directory scope where targets are created

# Compile flag overwrites
# if( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
#     set_source_files_properties( radiation_mcica_sw.F90 PROPERTIES COMPILE_OPTIONS "-vecabi=cmdtarget")
#     set_source_files_properties( radiation_cloud_generator.F90 PROPERTIES COMPILE_OPTIONS "-vecabi=cmdtarget")
# endif()

# For ecrad_ifs_driver_blocked we have to disable bounds checking
# (which is enabled in Debug builds) because sequence association
# and values initialized to -9999 for unused arguments is used in
# the call to RADIATION_SCHEME.
if( CMAKE_BUILD_TYPE MATCHES "Debug" )
    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        set_source_files_properties(
            "ecrad_ifs_driver_blocked.F90"
            PROPERTIES COMPILE_FLAGS "-fcheck=no-bounds"
        )
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        set_source_files_properties(
            "ecrad_ifs_driver_blocked.F90"
            PROPERTIES COMPILE_FLAGS "-check nobounds"
        )
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC")
        set_source_files_properties(
            "ecrad_ifs_driver_blocked.F90"
            PROPERTIES COMPILE_FLAGS "-Mnobounds"
        )
    endif()
endif()
