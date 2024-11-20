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
