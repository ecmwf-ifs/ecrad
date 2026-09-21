MODULE YOERRTO8

USE PARKIND1  ,ONLY : JPIM     ,JPRB,JPRD

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO8* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 8
!     BAND 8:  1080-1180 cm-1 (low (i.e.>~300mb) - H2O; high - O3)
!     ABozzo 201306 updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO8  = 16

REAL(KIND=JPRB) , DIMENSION(NO8) :: FRACREFAO
REAL(KIND=JPRB) , DIMENSION(NO8) :: FRACREFBO
REAL(KIND=JPRB) , DIMENSION(NO8) :: CFC12O
REAL(KIND=JPRB) , DIMENSION(NO8) :: CFC22ADJO

REAL(KIND=JPRB) :: KAO(NO8,5,13)
REAL(KIND=JPRB) :: KBO(NO8,5,13:59)
REAL(KIND=JPRD) :: KAO_D(NO8,5,13)
REAL(KIND=JPRD) :: KBO_D(NO8,5,13:59)
REAL(KIND=JPRB) :: SELFREFO(NO8,10)
REAL(KIND=JPRB) :: KAO_MCO2(NO8,19)
REAL(KIND=JPRB) :: KAO_MN2O(NO8,19)
REAL(KIND=JPRB) :: KAO_MO3(NO8,19)
REAL(KIND=JPRB) :: KBO_MCO2(NO8,19)
REAL(KIND=JPRB) :: KBO_MN2O(NO8,19)
REAL(KIND=JPRB) :: FORREFO(NO8,4)



!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSCO2A : REAL     
! ABSCO2B : REAL     
! ABSN2OA : REAL     
! ABSN2OB : REAL 
! CFC12   : REAL     
! CFC22ADJ: REAL     
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTO8
