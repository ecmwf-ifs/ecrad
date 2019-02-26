MODULE YOERRTAB

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

!    -------------------------------------------------------------------

REAL(KIND=JPRB) , DIMENSION(0:5000) :: TRANS
REAL(KIND=JPRB) :: BPADE

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! TRANS  :  REAL   : TABULATED REFERENCE FOR LW TRANSMISSION    
! BPADE  :  REAL   : INVERSE OF PADE APPROXIMATION CONSTANT  (= 1./0.278)     
!     -----------------------------------------------------------------
END MODULE YOERRTAB

