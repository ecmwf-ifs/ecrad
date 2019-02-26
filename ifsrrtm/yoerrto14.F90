MODULE YOERRTO14

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO14* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 14
!     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)
!     ABozzo 201306 updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO14 = 16

REAL(KIND=JPRB) , DIMENSION(NO14) :: FRACREFAO
REAL(KIND=JPRB) , DIMENSION(NO14) :: FRACREFBO

REAL(KIND=JPRB) :: KAO(5,13,NO14)
REAL(KIND=JPRB) :: KBO(5,13:59,NO14)
REAL(KIND=JPRD) :: KAO_D(5,13,NO14)
REAL(KIND=JPRD) :: KBO_D(5,13:59,NO14)
REAL(KIND=JPRB) :: SELFREFO(10,NO14)
REAL(KIND=JPRB) :: FORREFO(4,NO14)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTO14
