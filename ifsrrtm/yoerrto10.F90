MODULE YOERRTO10

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO14* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 10
!     BAND 10:  1390-1480 cm-1 (low - H2O; high - H2O)
!     ABozzo 201306 updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO10 = 16

REAL(KIND=JPRB) , DIMENSION(NO10) :: FRACREFAO
REAL(KIND=JPRB) , DIMENSION(NO10) :: FRACREFBO

REAL(KIND=JPRB) :: KAO(5,13,NO10)
REAL(KIND=JPRB) :: KBO(5,13:59,NO10)
REAL(KIND=JPRD) :: KAO_D(5,13,NO10)
REAL(KIND=JPRD) :: KBO_D(5,13:59,NO10)
REAL(KIND=JPRB) :: SELFREFO(10,NO10)
REAL(KIND=JPRB) :: FORREFO(4,NO10)

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
