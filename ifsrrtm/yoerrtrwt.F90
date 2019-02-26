MODULE YOERRTRWT

USE PARKIND1, ONLY : JPRB
!USE PARRRTM,  ONLY : JPGPT, JPG, JPBAND
USE PARRRTM,  ONLY : JPG, JPBAND, JPGMAX

IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

!    -------------------------------------------------------------------

REAL(KIND=JPRB) :: FREFA  (JPGMAX,13)
REAL(KIND=JPRB) :: FREFB  (JPGMAX,6)
REAL(KIND=JPRB) :: FREFADF(JPGMAX,13)
REAL(KIND=JPRB) :: FREFBDF(JPGMAX,6)
REAL(KIND=JPRB) :: RWGT   (JPG*JPBAND)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! FREFA  :  REAL   :
! FREFB  :  REAL   :
! FREFADF:  REAL   :
! FREFBDF:  REAL   :
! RWT    :  REAL   : 
!    -------------------------------------------------------------------
END MODULE YOERRTRWT
