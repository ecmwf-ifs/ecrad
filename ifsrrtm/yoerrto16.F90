MODULE YOERRTO16

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO16* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 16
!     BAND 16:  2600-3000 cm-1 (low - H2O,CH4; high - nothing)
!     ABozzo 201306 updated to rrtmg v4.85
!     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO16 = 16

REAL(KIND=JPRB) :: FRACREFAO(NO16,9)
REAL(KIND=JPRB) , DIMENSION(NO16) :: FRACREFBO

REAL(KIND=JPRB) :: KAO(9,5,13,NO16)
REAL(KIND=JPRB) :: KBO(5,13:59,NO16)
REAL(KIND=JPRD) :: KAO_D(9,5,13,NO16)
REAL(KIND=JPRD) :: KBO_D(5,13:59,NO16)
REAL(KIND=JPRB) :: SELFREFO(10,NO16)
REAL(KIND=JPRB) :: FORREFO(4,NO16)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTO16
