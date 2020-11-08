/**
 * (C) Copyright 2014- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

INTERFACE
SUBROUTINE EC_MPI_FINALIZE(KERROR,LDCALLFINITO,LDMEMINFO,CALLER)
USE PARKIND1, ONLY : JPIM
INTEGER(KIND=JPIM), INTENT(OUT) :: KERROR
LOGICAL, INTENT(IN) :: LDCALLFINITO
LOGICAL, INTENT(IN) :: LDMEMINFO
CHARACTER(LEN=*), INTENT(IN) :: CALLER
END SUBROUTINE EC_MPI_FINALIZE
END INTERFACE
