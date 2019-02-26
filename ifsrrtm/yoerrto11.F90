MODULE YOERRTO11

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO11* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 11
!     BAND 11:  1480-1800 cm-1 (low - H2O; high - H2O)
!     ABozzo 201306 updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO11 = 16

REAL(KIND=JPRB) , DIMENSION(NO11) :: FRACREFAO
REAL(KIND=JPRB) , DIMENSION(NO11) :: FRACREFBO

REAL(KIND=JPRB) :: KAO(5,13,NO11)
REAL(KIND=JPRB) :: KBO(5,13:59,NO11)
REAL(KIND=JPRD) :: KAO_D(5,13,NO11)
REAL(KIND=JPRD) :: KBO_D(5,13:59,NO11)
REAL(KIND=JPRB) :: KAO_MO2(19,NO11)
REAL(KIND=JPRB) :: KBO_MO2(19,NO11)
REAL(KIND=JPRB) :: SELFREFO(10,NO11)
REAL(KIND=JPRB) :: FORREFO(4,NO11)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTO11
