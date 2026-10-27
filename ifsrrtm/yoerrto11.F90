MODULE YOERRTO11

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO11* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 11
!     BAND 11:  1480-1800 cm-1 (low - H2O; high - H2O)
!     ABozzo 201306 updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO11 = 16

REAL(KIND=JPRB) , DIMENSION(NO11) :: FRACREFAO
REAL(KIND=JPRB) , DIMENSION(NO11) :: FRACREFBO

REAL(KIND=JPRB) :: KAO(NO11,5,13)
REAL(KIND=JPRB) :: KBO(NO11,5,13:59)
REAL(KIND=JPRD) :: KAO_D(NO11,5,13)
REAL(KIND=JPRD) :: KBO_D(NO11,5,13:59)
REAL(KIND=JPRB) :: KAO_MO2(NO11,19)
REAL(KIND=JPRB) :: KBO_MO2(NO11,19)
REAL(KIND=JPRB) :: SELFREFO(NO11,10)
REAL(KIND=JPRB) :: FORREFO(NO11,4)

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
