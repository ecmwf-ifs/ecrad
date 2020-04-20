MODULE MPL_SETDFLT_COMM_MOD

!**** MPL_SETDFLT_COMM Set new default communicator

!     Purpose.
!     --------
!     Set new communicator as default, and return old communicator

!**   Interface.
!     ----------
!        CALL MPL_SETDFLT_COMM(KCOMM,KCOMM_OLD)

!        Input required arguments :
!        -------------------------
!           KCOMM    -  New communicator

!        Input optional arguments :
!        -------------------------

!        Output required arguments :
!        -------------------------
!           KCOMM_OLD    -  Old communicator

!        Output optional arguments :
!        -------------------------

!     Author.
!     -------
!        J.Hague        

!     Modifications.
!     --------------
!        Original: 2003-16-07

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE MPL_MPIF
USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD

IMPLICIT NONE

PRIVATE
PUBLIC MPL_SETDFLT_COMM

CONTAINS 


SUBROUTINE MPL_SETDFLT_COMM(KCOMM,KCOMM_OLD)


#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_COMM_SIZE => MPI_COMM_SIZE8
#endif


INTEGER(KIND=JPIM),INTENT(IN)  :: KCOMM
INTEGER(KIND=JPIM),INTENT(OUT) :: KCOMM_OLD

INTEGER(KIND=JPIM) :: IER
INTEGER(KIND=JPIM) :: ITID
ITID = OML_MY_THREAD()

KCOMM_OLD=MPL_COMM_OML(ITID)
MPL_COMM_OML(ITID)=KCOMM

CALL MPI_COMM_SIZE(KCOMM,MPL_NUMPROC,IER)

RETURN
END SUBROUTINE MPL_SETDFLT_COMM

END MODULE MPL_SETDFLT_COMM_MOD
