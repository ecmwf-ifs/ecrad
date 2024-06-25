# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

set(
    utilities_SOURCES
        easy_netcdf_read_mpi.F90
        easy_netcdf.F90
        radiation_io.F90
        radiation_random_numbers_mix.F90
)

ecbuild_add_library(
    TARGET ecrad_utilities.${PREC}
    TYPE OBJECT
    SOURCES ${utilities_SOURCES}
    PUBLIC_INCLUDES
        "$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/radiation.${PREC}>"
    PRIVATE_LIBS
        ecrad_base.${PREC}
        ecrad_ifsaux.${PREC}
    PUBLIC_LIBS
        NetCDF::NetCDF_Fortran
)

install(
    FILES
        "${CMAKE_Fortran_MODULE_DIRECTORY}/easy_netcdf_read_mpi.mod"
        "${CMAKE_Fortran_MODULE_DIRECTORY}/easy_netcdf.mod"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/radiation.${PREC}"
)
