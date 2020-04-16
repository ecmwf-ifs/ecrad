INTERFACE
SUBROUTINE USER_CLOCK(PELAPSED_TIME,PELAPSED_TIME_SINCE,PVECTOR_CP,PTOTAL_CP)

!**** *USER_CLOCK* - interface to system dependent timer routines

!     Purpose.
!     --------
!        Returns elapsed and CP from the start of execution.
!        Elapsed time is made relative to the first call to USER_CLOCK.

!**   Interface.
!     ----------
!        ZTIME=USER_CLOCK(PELAPSED_TIME,PELAPSED_TIME_SINCE,
!                         PVECTOR_CP,PTOTAL_CP)

!        Explicit arguments: (All are optional arguments)
!                           PELAPSED_TIME=wall clock time (seconds)
!                           PELAPSED_TIME_SINCE=wall clock time (seconds)
!                             change from input value of this parameter
!                           PVECTOR_CP=CP vector time  (seconds)
!                           PTOTAL_CP=total CP time   (seconds)

!     Author.
!     -------
!        D.Dent      *ECMWF*

!     External References:
!     -------------------

!        TIMEF,CPTIME

!     Modifications.
!     --------------
!        Original  : 97-09-25
!      F. Vana  05-Mar-2015  Support for single precision
!     ----------------------------------------------------------


USE PARKIND1  ,ONLY : JPRD, JPIM

IMPLICIT NONE

REAL(KIND=JPRD),INTENT(OUT) :: PELAPSED_TIME,PVECTOR_CP,PTOTAL_CP
REAL(KIND=JPRD),INTENT(INOUT) :: PELAPSED_TIME_SINCE
OPTIONAL            PELAPSED_TIME,PELAPSED_TIME_SINCE
OPTIONAL            PVECTOR_CP,PTOTAL_CP
REAL(KIND=JPRD)      :: ZVECTOR_CP,ZTOTAL_CP,ZWALL
REAL(KIND=JPRD),EXTERNAL :: TIMEF


END SUBROUTINE USER_CLOCK
END INTERFACE
