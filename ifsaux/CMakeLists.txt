# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

set( ifsaux_SOURCES
        yomcst.F90
        yomdyncore_ecrad.F90
        yomlun_ifsaux.F90
        yomlun_ecrad.F90
        yommp0_ifsaux.F90
        yomtag.F90
)

if( NOT HAVE_FIAT )
    list( APPEND ifsaux_SOURCES parkind1.F90 mpl_module.F90 abor1.F90 )
endif()

ecbuild_add_library(
    TARGET ecrad_ifsaux.${PREC}
    TYPE OBJECT
    SOURCES ${ifsaux_SOURCES}
    PUBLIC_INCLUDES
        "$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>"
    PRIVATE_LIBS ecrad_base.${PREC}
    PRIVATE_DEFINITIONS ${GPU_OFFLOAD}GPU
)
install( TARGETS ecrad_ifsaux.${PREC} EXPORT ecrad-targets )
