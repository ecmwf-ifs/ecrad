SUBROUTINE GETMEMVALS(N, KEY, KVAL)
USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPIB
IMPLICIT NONE
INTEGER(KIND=JPIM), INTENT(IN) :: N, KEY(N)
INTEGER(KIND=JPIB), INTENT(OUT):: KVAL(N)
!--------------------------------- key ----------------------------------------------
INTEGER(KIND=JPIB), EXTERNAL :: GETHWM    !  1  High Water Mark for HEAP-alloc
INTEGER(KIND=JPIB), EXTERNAL :: GETMAXRSS !  2  Maximum resident memory so far
INTEGER(KIND=JPIB), EXTERNAL :: GETCURHEAP!  3  Instantaneous allocation from ALLOCATE/malloc
INTEGER(KIND=JPIB), EXTERNAL :: GETSTK    !  4  Instantaneous stack usage
INTEGER(KIND=JPIB), EXTERNAL :: GETMAXSTK !  5  Maximum stack usage so far
INTEGER(KIND=JPIB), EXTERNAL :: GETPAG    !  6  I/O caused by paging
! -- add more as required (all 64-bit integers upon return, though) --

INTEGER(KIND=JPIM) J

DO J=1,N
  IF (KEY(J) == 1) THEN
    KVAL(J) = GETHWM()
  ELSE IF (KEY(J) == 2) THEN
    KVAL(J) = GETMAXRSS()
  ELSE IF (KEY(J) == 3) THEN
    KVAL(J) = GETCURHEAP()
  ELSE IF (KEY(J) == 4) THEN
    KVAL(J) = GETSTK()
  ELSE IF (KEY(J) == 5) THEN
    KVAL(J) = GETMAXSTK()
  ELSE IF (KEY(J) == 6) THEN
    KVAL(J) = GETPAG()
  ENDIF
ENDDO

END SUBROUTINE GETMEMVALS
