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
SUBROUTINE DR_HOOK_UTIL_MULTI(LDHOOK,CDNAME,KCASE,PKEY,KPKEY,CDFILENAME,KSIZEINFO)
USE PARKIND1  ,ONLY : JPIM     ,JPRB
IMPLICIT NONE
LOGICAL,INTENT(IN)          :: LDHOOK
CHARACTER(LEN=*),INTENT(IN) :: CDNAME,CDFILENAME
INTEGER(KIND=JPIM),INTENT(IN) :: KPKEY, KCASE,KSIZEINFO
REAL(KIND=JPRB),INTENT(INOUT) :: PKEY(KPKEY)
END SUBROUTINE DR_HOOK_UTIL_MULTI
END INTERFACE
