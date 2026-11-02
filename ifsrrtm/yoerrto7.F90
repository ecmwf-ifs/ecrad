MODULE YOERRTO7

USE PARKIND1  ,ONLY : JPIM     ,JPRB,JPRD

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO7* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 7
!     BAND 7:  980-1080 cm-1 (low - H2O,O3; high - O3)
!     ABozzo updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO7  = 16

REAL(KIND=JPRB) :: FRACREFAO(NO7,9)

REAL(KIND=JPRB) , DIMENSION(NO7) :: FRACREFBO
REAL(KIND=JPRB) :: KAO(NO7,9,5,13)
REAL(KIND=JPRB) :: KBO(NO7,5,13:59)
REAL(KIND=JPRD) :: KAO_D(NO7,9,5,13)
REAL(KIND=JPRD) :: KBO_D(NO7,5,13:59)
REAL(KIND=JPRB) :: SELFREFO(NO7,10)
REAL(KIND=JPRB) :: KAO_MCO2(NO7,9,19)
REAL(KIND=JPRB) :: KBO_MCO2(NO7,19)
REAL(KIND=JPRB) :: FORREFO(NO7,4)



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
END MODULE YOERRTO7
