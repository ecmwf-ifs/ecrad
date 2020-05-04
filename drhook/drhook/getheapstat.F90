SUBROUTINE GETHEAPSTAT(KOUT, CDLABEL)

USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPIB

USE MPL_MODULE

#ifdef NAG
USE F90_UNIX_ENV, ONLY: GETENV
#endif

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KOUT
CHARACTER(LEN=*), INTENT(IN) :: CDLABEL
INTEGER(KIND=JPIM) :: I, IMYPROC, INPROC, IRET, IOFFSET, II
INTEGER(KIND=JPIM), PARAMETER :: JP_NPROFILE = 9 ! pls. consult ifsaux/utilities/getcurheap.c
INTEGER(KIND=JPIM), PARAMETER :: ISIZE = JP_NPROFILE+1
INTEGER(KIND=JPIB) ILIMIT(ISIZE)
INTEGER(KIND=JPIB) ICNT(ISIZE)
REAL(KIND=JPRB), ALLOCATABLE :: ZSEND(:), ZRECV(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: ICOUNTS(:)
CHARACTER(LEN=1) CLENV
CHARACTER(LEN=80) CLTEXT(0:4)

CALL GET_ENVIRONMENT_VARIABLE("EC_PROFILE_HEAP", CLENV) ! turn OFF by export EC_PROFILE_HEAP=0

IF (KOUT >= 0 .AND. CLENV /= '0') THEN
  IMYPROC = MPL_MYRANK()
  INPROC  = MPL_NPROC()

  DO I=1,ISIZE
    ILIMIT(I) = I ! power of 10's ; pls. consult ifsaux/utilities/getcurheap.c
  ENDDO

  ALLOCATE(ZSEND(ISIZE))
  ALLOCATE(ZRECV(ISIZE * INPROC))
  ALLOCATE(ICOUNTS(INPROC))

  CLTEXT(0) = "free()/DEALLOCATE -hits per byte range"
  CLTEXT(1) = "malloc()/ALLOCATE -hits per byte range"
  CLTEXT(2) = "Outstanding malloc()/ALLOCATE -hits per byte range"
  CLTEXT(3) = "Outstanding amount of malloc()/ALLOCATE -bytes per byte range"
  CLTEXT(4) = "Average amount of outstanding malloc()/ALLOCATE -bytes per byte range"

  DO II=0,4
    ICNT(:) = 0
    CALL PROFILE_HEAP_GET(ICNT, ISIZE, II, IRET)

    ZSEND(:) = 0
    DO I=1,IRET
      ZSEND(I) = ICNT(I)
    ENDDO
    ZRECV(:) = -1

    ICOUNTS(:) = ISIZE
    CALL MPL_GATHERV(ZSEND(:), KROOT=1, KRECVCOUNTS=ICOUNTS(:), &
                    &PRECVBUF=ZRECV, CDSTRING='GETHEAPSTAT:')

    IF (IMYPROC == 1) THEN
!     Not more than 132 columns, please :-)
      WRITE(KOUT,9000) TRIM(CLTEXT(II)),TRIM(CDLABEL), "Node", &
                     & (ILIMIT(I),I=1,MIN(JP_NPROFILE,9)), "Larger"
9000  FORMAT(/,"Heap Utilization Profile (",A,"): ",A,&
            &/,126("="),&
            &//,(A4,2X,9(:,2X,4X,"< 10^",I1),:,2X,A10))
      WRITE(KOUT,9001)
9001  FORMAT(4("="),2X,10(2X,10("="))/)
      IOFFSET = 0
      DO I=1,INPROC
        ICNT(:) = ZRECV(IOFFSET+1:IOFFSET+ISIZE)
        WRITE(KOUT,'(i4,2x,(10(:,2x,i10)))') I,ICNT(:)
        IOFFSET = IOFFSET + ISIZE
      ENDDO
    ENDIF
  ENDDO

  IF (IMYPROC == 1) THEN
    WRITE(KOUT,'(/,a,/)') 'End of Heap Utilization Profile'
  ENDIF

  DEALLOCATE(ZSEND)
  DEALLOCATE(ZRECV)
  DEALLOCATE(ICOUNTS)
ENDIF
END SUBROUTINE GETHEAPSTAT
