MODULE YOMMPI


USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    Contains identifiers used by MPI (Message Passing Interface)


INTEGER(KIND=JPIM), PARAMETER :: MINTET = 1
INTEGER(KIND=JPIM), PARAMETER :: MREALT = 2
INTEGER(KIND=JPIM), PARAMETER :: MLOGIT = 3
INTEGER(KIND=JPIM), PARAMETER :: MCHART = 4

END MODULE YOMMPI
