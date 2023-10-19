! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

SUBROUTINE DR_HOOK(CDNAME,KSWITCH,PKEY)
USE PARKIND1, ONLY : JPIM, JPRD, JPRB
USE YOMHOOK, ONLY: LHOOK
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: CDNAME
INTEGER(KIND=JPIM),        INTENT(IN) :: KSWITCH
REAL(KIND=JPRD),        INTENT(INOUT) :: PKEY
REAL(KIND=JPRB) :: ZKEY
ZKEY = TRANSFER(PKEY,ZKEY)
CALL DR_HOOK_UTIL(LHOOK,CDNAME,KSWITCH,ZKEY,'',0_JPIM)
PKEY = TRANSFER(ZKEY,PKEY)
END SUBROUTINE DR_HOOK
