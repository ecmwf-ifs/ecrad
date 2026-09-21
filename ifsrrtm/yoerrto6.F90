MODULE YOERRTO6

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO6* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 6
!     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)
!     ABozzo 201306 update to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO6  = 16

REAL(KIND=JPRB) , DIMENSION(NO6) :: FRACREFAO

REAL(KIND=JPRB) , DIMENSION(NO6) :: CFC11ADJO
REAL(KIND=JPRB) , DIMENSION(NO6) :: CFC12O


REAL(KIND=JPRB) :: KAO(NO6,5,13)
REAL(KIND=JPRD) :: KAO_D(NO6,5,13)
REAL(KIND=JPRB) :: SELFREFO(NO6,10)
REAL(KIND=JPRB) :: KAO_MCO2(NO6,19)
REAL(KIND=JPRB) :: FORREFO(NO6,4)



!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTO6
