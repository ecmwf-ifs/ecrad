MODULE YOERRTO13

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO13* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 13
!     BAND 13:  2080-2250 cm-1 (low - H2O,N2O; high - nothing)
!     ABozzo 201306 updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO13 = 16

REAL(KIND=JPRB) :: FRACREFAO(NO13,9)
REAL(KIND=JPRB) , DIMENSION(NO13) :: FRACREFBO

REAL(KIND=JPRB) :: KAO(9,5,13,NO13)
REAL(KIND=JPRD) :: KAO_D(9,5,13,NO13)
REAL(KIND=JPRB) :: SELFREFO(10,NO13)
REAL(KIND=JPRB) :: KAO_MCO2(9,19,NO13)
REAL(KIND=JPRB) :: KAO_MCO(9,19,NO13)
REAL(KIND=JPRB) :: KBO_MO3(19,NO13)
REAL(KIND=JPRB) :: FORREFO(4,NO13)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL
!     -----------------------------------------------------------------
END MODULE YOERRTO13
