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
)

ecbuild_add_library(
    TARGET ifs.${PREC}
    TYPE OBJECT
    SOURCES ${ifs_SOURCES}
    PRIVATE_LIBS
        ecrad.${PREC}
        ecrad_base.${PREC}
    PUBLIC_LIBS
        $<${HAVE_NVTX}:${NVTX_TARGET}>
)
