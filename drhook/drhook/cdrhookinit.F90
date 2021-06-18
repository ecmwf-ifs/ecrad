! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

SUBROUTINE CDRHOOKINIT(KRET)
!-- Makes sure Dr.Hook gets properly initialized from C-main program, too
USE PARKIND1  ,ONLY  : JPIM, JPRB, JPRD
USE YOMHOOK   ,ONLY  : LHOOK, DR_HOOK
IMPLICIT NONE
INTEGER(KIND=JPIM), INTENT(OUT) :: KRET
REAL(KIND=JPRD) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('CDRHOOKINIT',0,ZHOOK_HANDLE)
IF (LHOOK) THEN
  KRET = 1 ! Dr.Hook is ON
ELSE
  KRET = 0 ! Dr.Hook is OFF
ENDIF
IF (LHOOK) CALL DR_HOOK('CDRHOOKINIT',1,ZHOOK_HANDLE)
END SUBROUTINE CDRHOOKINIT
