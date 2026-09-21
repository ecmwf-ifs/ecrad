MODULE YOERRTO10

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO14* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 10
!     BAND 10:  1390-1480 cm-1 (low - H2O; high - H2O)
!     ABozzo 201306 updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO10 = 16

REAL(KIND=JPRB) , DIMENSION(NO10) :: FRACREFAO
REAL(KIND=JPRB) , DIMENSION(NO10) :: FRACREFBO

REAL(KIND=JPRB) :: KAO(NO10,5,13)
REAL(KIND=JPRB) :: KBO(NO10,5,13:59)
REAL(KIND=JPRD) :: KAO_D(NO10,5,13)
REAL(KIND=JPRD) :: KBO_D(NO10,5,13:59)
REAL(KIND=JPRB) :: SELFREFO(NO10,10)
REAL(KIND=JPRB) :: FORREFO(NO10,4)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTO10
