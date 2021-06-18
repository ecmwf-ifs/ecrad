! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

SUBROUTINE ABOR1(CDTEXT)

USE PARKIND1,      ONLY : JPIM, JPRB
USE YOMLUN_IFSAUX, ONLY : NULOUT, NULERR

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: CDTEXT

IF (NULOUT >= 0) WRITE(NULOUT,'(1X,A)') CDTEXT
IF (NULERR >= 0) WRITE(NULERR,'(1X,A,A)') 'ABORT! ', CDTEXT

IF (NULOUT >= 0) THEN
  ! FLUSH not understood by NAG compiler
  !CALL FLUSH(NULOUT)
  IF (NULOUT /= 0 .and. NULOUT /= 6) CLOSE(NULOUT)
ENDIF

#ifdef __PGI
      stop 1
#else
      error stop 1
#endif

END SUBROUTINE ABOR1
