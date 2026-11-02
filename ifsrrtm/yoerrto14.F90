MODULE YOERRTO14

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO14* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 14
!     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)
!     ABozzo 201306 updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO14 = 16

REAL(KIND=JPRB) , DIMENSION(NO14) :: FRACREFAO
REAL(KIND=JPRB) , DIMENSION(NO14) :: FRACREFBO

REAL(KIND=JPRB) :: KAO(NO14,5,13)
REAL(KIND=JPRB) :: KBO(NO14,5,13:59)
REAL(KIND=JPRD) :: KAO_D(NO14,5,13)
REAL(KIND=JPRD) :: KBO_D(NO14,5,13:59)
REAL(KIND=JPRB) :: SELFREFO(NO14,10)
REAL(KIND=JPRB) :: FORREFO(NO14,4)

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
