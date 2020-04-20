#ifdef RS6K

!-- Thanks to John Hague (IBM) and Deborah Salmond (ECMWF) -- 2002
!   real --> real*8 by Sami Saarinen, ECMWF, 16-Jan-2004

      SUBROUTINE JDSORT(N, ISTRIDE, IPTR, C, INDEX)
!
!**C* *JDSORT* - interface to IBM essl dsorts
!
!**C  INTERFACE.
!     ----------
!
!          *CALL* *JDSORT(N,ISTRIDE, IPTR, C, INDEX)
!           N = number of elements to be sorted
!           ISTRIDE = stride
!           IPTR = start address
!           C = array to be sorted
!           INDEX = index into array
!
      USE PARKIND1, ONLY : JPRB

      IMPLICIT NONE
      INTEGER N,ISTRIDE, IPTR
      REAL(KIND=JPRB) :: C(*)
      INTEGER INDEX(N)

! local variables

      INTEGER I,J,IXX(N),IX(N)
      REAL(KIND=JPRB) :: CC(N)
      REAL(KIND=JPRB) :: W(N/2+1)

      DO I=1,N
        J=(I-1)*ISTRIDE+1
        CC(I)=C(IPTR+(INDEX(I)-1)*ISTRIDE)
        IXX(I)=INDEX(I)
      ENDDO
      CALL DSORTS(CC,1,N,IX,W,SIZE(W))
      DO I=1,N
        INDEX(I)=IXX(IX(I))
      ENDDO
      RETURN
      END SUBROUTINE JDSORT


      SUBROUTINE JISORT(N, ISTRIDE, IPTR, C, INDEX)
!
!**** *JISORT* - interface to IBM essl isorts
!
!**   INTERFACE.
!     ----------
!
!          *CALL* *JISORT(N,ISTRIDE, IPTR, C, INDEX)
!           N = number of elements to be sorted
!           ISTRIDE = stride
!           IPTR = start address
!           C = array to be sorted
!           INDEX = index into array
!
      IMPLICIT NONE
      INTEGER N,ISTRIDE, IPTR
      INTEGER C(*),INDEX(N)

! local variables

      INTEGER CC(N)
      INTEGER I,J,IXX(N),IX(N)
      INTEGER W(N/2+1)

      DO I=1,N
        J=(I-1)*ISTRIDE+1
        CC(I)=C(IPTR+(INDEX(I)-1)*ISTRIDE)
        IXX(I)=INDEX(I)
      ENDDO
      CALL ISORTS(CC,1,N,IX,W,SIZE(W))
      DO I=1,N
        INDEX(I)=IXX(IX(I))
      ENDDO
      RETURN
      END SUBROUTINE JISORT
#else
!-- other machines --
      SUBROUTINE JSORT_DUMMY()
      END SUBROUTINE JSORT_DUMMY
#endif

