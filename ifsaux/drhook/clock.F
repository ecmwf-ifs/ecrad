      FUNCTION CLOCK()
!
!     clock() returns the amount of CPU time (in microseconds) used since
!     the first call to clock().  The time reported is the sum of the user
!     and system times of the calling process and its terminated child
!     processes for which it has executed wait() or system() (see wait(2)
!     and system(3S)).  To determine the time in seconds, the value returned
!     by clock() should be divided by the value of the macro CLOCKS_PER_SEC.
!
!     INTRODUCE A PROPER CLCOK ROUTINE HERE
!
!      IMPLICIT LOGICAL (L)
      USE PARKIND1, ONLY : JPIM, JPRB
!
      IMPLICIT NONE
!
      INTEGER(KIND=JPIM) :: CLOCK
      INTEGER(KIND=JPIM) :: IC, ICR, ICM
      CALL SYSTEM_CLOCK(IC,ICR,ICM)
      IF(ICM.EQ.0) THEN
        PRINT *,' no system clock'
        CLOCK=-99
      ELSE
#ifdef VPP
        CLOCK=IC
#else
        CLOCK=INT(1000.0_JPRB*REAL(IC,KIND=JPRB)/ICR)
#endif
      ENDIF
      RETURN
      ENDFUNCTION CLOCK
