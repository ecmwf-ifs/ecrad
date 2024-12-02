# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

set( CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )

ecbuild_add_library(
    TARGET drhook_dummy
    SOURCES yomhook_dummy.F90
    PUBLIC_INCLUDES "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
)
