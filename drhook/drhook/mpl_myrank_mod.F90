MODULE MPL_MYRANK_MOD

!**** MPL_MYRANK - Find rank

!     Purpose.
!     --------
!     Returns the rank of the calling process 
!     in the currently active communicator

!**   Interface.
!     ----------
!        IRANK=MPL_MYRANK(KCOMM)

!        Input optional arguments :
!        -------------------------
!           KCOMM    -  Communicator number if different from MPI_COMM_WORLD 

!        Input required arguments :
!        -------------------------

!        Input optional arguments :
!        -------------------------
!           none

!        Output required arguments :
!        -------------------------
!           none

!        Output optional arguments :
!        -------------------------
!           none
!     Author.
!     -------
!        D.Dent, M.Hamrud     ECMWF

!     Modifications.
!     --------------
!        Original: 2000-09-01

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD

IMPLICIT NONE
PRIVATE
PUBLIC MPL_MYRANK

CONTAINS 

FUNCTION MPL_MYRANK(KCOMM)

#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_COMM_RANK => MPI_COMM_RANK8
#endif

INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KCOMM
INTEGER(KIND=JPIM) :: MPL_MYRANK

INTEGER(KIND=JPIM) IRANK,IERROR,ICOMM
LOGICAL :: LLABORT=.TRUE.

IF(MPL_NUMPROC < 1) CALL MPL_MESSAGE( &
  & CDMESSAGE='MPL_MYRANK: MPL NOT INITIALISED ',LDABORT=LLABORT) 

IF(PRESENT(KCOMM)) THEN
  CALL MPI_COMM_RANK(KCOMM, IRANK, IERROR)
  IF(IERROR /= 0 ) CALL MPL_MESSAGE(IERROR,&
   &'MPL_MYRANK: ERROR FROM MPI_COMM_RANK')
  MPL_MYRANK=IRANK+1
ELSE
  MPL_MYRANK=MPL_RANK
ENDIF

END FUNCTION MPL_MYRANK

END MODULE MPL_MYRANK_MOD
