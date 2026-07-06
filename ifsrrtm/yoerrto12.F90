MODULE YOERRTO12

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO12* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 12
!     BAND 12:  1800-2080 cm-1 (low - H2O,CO2; high - nothing)
!     ABozzo 201306 updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO12 = 16

REAL(KIND=JPRB) :: FRACREFAO(NO12,9)
REAL(KIND=JPRB) :: KAO(NO12,9,5,13)
REAL(KIND=JPRD) :: KAO_D(NO12,9,5,13)
REAL(KIND=JPRB) :: SELFREFO(NO12,10)
REAL(KIND=JPRB) :: FORREFO(NO12,4)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL
!     -----------------------------------------------------------------
END MODULE YOERRTO12
