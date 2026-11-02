# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

ecbuild_add_test(
    TARGET ecrad_${PREC}_cloud_cover_ranges
    SOURCES cloud_cover_ranges.F90
    LIBS ecrad.${PREC} ecrad_base.${PREC}
)

# Keep the last precision from the top-level build loop out of NVHPC's
# module search path: its module output directory takes precedence over -I.
set_target_properties( ecrad_${PREC}_cloud_cover_ranges PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module" )
file( MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module" )
