! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

SUBROUTINE ABOR1FL(CDFILE, KLINENUM, CDTEXT)
USE PARKIND1  ,ONLY : JPIM
USE YOMLUN_IFSAUX, ONLY : NULOUT,NULERR
#ifdef NAG
USE F90_UNIX_IO, ONLY: FLUSH
#endif
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: CDFILE,CDTEXT
INTEGER(KIND=JPIM), INTENT(IN) :: KLINENUM
IF (LEN(CDFILE) > 0 .AND. KLINENUM > 0 .AND. NULERR >= 0) THEN
 1000 FORMAT(1X,A,A,":",I6.6)
  WRITE(NULERR,1000) 'ABOR1FL HAS BEEN CALLED AT ',CDFILE,KLINENUM
  CALL FLUSH(NULERR)
ENDIF
CALL ABOR1(CDTEXT)
END SUBROUTINE ABOR1FL
