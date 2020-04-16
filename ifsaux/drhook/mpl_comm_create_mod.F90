MODULE MPL_COMM_CREATE_MOD

!**** MPL_COMM_CREATE Create a new communicator

!     Purpose.
!     --------
!     Create a new communicator and set as default

!**   Interface.
!     ----------
!        CALL MPL_COMM_CREATE

!        Input required arguments :
!        -------------------------

!        Input optional arguments :
!        -------------------------

!        Output required arguments :
!        -------------------------

!        Output optional arguments :
!        -------------------------
!           KERROR   -  return error code.     If not supplied, 
!                       MPL_COMM_CREATE aborts when an error is detected.
!     Author.
!     -------
!        D.Dent, M.Hamrud     ECMWF

!     Modifications.
!     --------------
!        Original: 2000-09-01

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE MPL_MPIF
USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD

IMPLICIT NONE

PRIVATE
PUBLIC MPL_COMM_CREATE

CONTAINS 

SUBROUTINE MPL_COMM_CREATE(KERROR)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KERROR
INTEGER(KIND=JPIM) :: ITID
ITID = OML_MY_THREAD()
!   this line to be replaced
MPL_COMM_OML(ITID)=MPL_COMM
KERROR=0

RETURN
END SUBROUTINE MPL_COMM_CREATE

END MODULE MPL_COMM_CREATE_MOD
