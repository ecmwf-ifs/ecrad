MODULE YOERRTO15

USE PARKIND1  ,ONLY : JPIM     ,JPRB,JPRD

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO15* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 15
!     BAND 15:  2380-2600 cm-1 (low - N2O,CO2; high - nothing)
!     ABozzo 2001306 updated to rrtmg v4.85
!     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
!                              (high - nothing)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO15 = 16

REAL(KIND=JPRB) :: FRACREFAO(NO15,9)

REAL(KIND=JPRB) :: KAO(9,5,13,NO15)
REAL(KIND=JPRD) :: KAO_D(9,5,13,NO15)
REAL(KIND=JPRB) :: SELFREFO(10,NO15)
REAL(KIND=JPRB) :: FORREFO(4,NO15)
REAL(KIND=JPRB) :: KAO_MN2(9,19,NO15)
!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL 
!     -----------------------------------------------------------------
END MODULE YOERRTO15
